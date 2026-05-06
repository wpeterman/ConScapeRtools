#' Inspect ConScape slot names used by ConScapeRtools
#'
#' @description
#' Returns a compact crosswalk between the R-facing argument names in this
#' package and the Julia-side ConScape slots/functions they populate.
#'
#' @details
#' ConScapeRtools keeps the original argument names (`hab_target`, `hab_src`,
#' `mov_prob`, and `exp_d`) for backwards compatibility, but new workflows
#' should prefer the clearer ConScape-aligned names shown in the `run_conscape`
#' column. In ConScape terminology, cell habitat rasters are source and target
#' qualities, while the movement raster is an affinity/permeability surface
#' that is converted into a graph adjacency matrix. Costs are then either
#' derived from affinities with a named transformation such as `"minuslog"` or
#' supplied by a custom Julia transformation.
#'
#' @return A data frame describing each API slot.
#' @export
conscape_api_slots <- function() {
  data.frame(
    slot = c(
      "target_qualities",
      "source_qualities",
      "affinities",
      "costs",
      "theta",
      "distance_scale",
      "landmark",
      "connectivity_function",
      "metrics",
      "sensitivity"
    ),
    run_conscape = c(
      "target_qualities (alias: hab_target)",
      "source_qualities (alias: hab_src)",
      "affinities (alias: mov_prob)",
      "cost_function",
      "theta",
      "distance_scale (alias: exp_d)",
      "landmark",
      "connectivity_function",
      "metrics",
      "sensitivity"
    ),
    conscape_prep = c(
      "r_target",
      "r_src",
      "r_mov",
      "cost_function in run_conscape()",
      "theta in run_conscape()",
      "distance_scale in run_conscape()",
      "landmark",
      "connectivity_function in run_conscape()",
      "metrics in run_conscape()",
      "sensitivity in run_conscape()"
    ),
    julia = c(
      "ConScape.Grid(target_qualities = ...)",
      "ConScape.Grid(source_qualities = ...)",
      "ConScape.Grid(affinities = ConScape.graph_matrix_from_raster(...))",
      "ConScape.Grid(costs = ...)",
      "ConScape.GridRSP(theta = ...)",
      "distance_transformation = x -> exp(-x / distance_scale)",
      "ConScape.coarse_graining(g, landmark)",
      "ConScape.expected_cost / least_cost_distance / free_energy_distance / survival_probability / power_mean_proximity",
      "ConScape.connected_habitat / betweenness_kweighted / betweenness_qweighted / criticality",
      "ConScape.sensitivity or ConScape.sensitivity_simulation"
    ),
    stringsAsFactors = FALSE
  )
}

#' Configure ConScape sensitivity outputs
#'
#' @description
#' Creates a sensitivity specification for [run_conscape()]. Analytical
#' sensitivity is available when the installed Julia ConScape package provides
#' `ConScape.sensitivity`, currently on the ConScape `sensitivity` branch or a
#' future release containing that function.
#'
#' @details
#' Sensitivity analysis asks how a landscape-level connectivity summary changes
#' when a local component of the landscape is perturbed. The `wrt` argument
#' controls which component is perturbed:
#'
#' * `"Q"` evaluates sensitivity to habitat quality.
#' * `"A"` evaluates sensitivity to affinities/permeability.
#' * `"C"` evaluates sensitivity to movement costs.
#' * `"A&C=f(A)"` evaluates the linked case where changing affinity also changes
#'   cost through the chosen cost transformation.
#' * `"C&A=f(C)"` evaluates the reverse linked case.
#'
#' The default `unitless = TRUE` returns elasticities: sensitivities scaled by
#' the value of the perturbed quantity. These are usually easier to compare
#' across habitat-quality and movement surfaces because they describe response
#' to proportional change. Set `unitless = FALSE` for raw unit-change
#' sensitivities.
#'
#' For tiled workflows, sensitivity is computed independently within each tile
#' and then mosaicked by [run_conscape()]. This is most defensible with
#' `landscape_measure = "sum"`, `landmark = 1L`, identical source and target
#' qualities, and a conservative tile overlap. The current ConScape
#' implementation assumes source and target qualities are equal; therefore
#' [run_conscape()] requires `landmark = 1L` by default when sensitivity is
#' requested so that target qualities are not altered by coarse graining.
#' Tile-wise `landscape_measure = "eigenanalysis"` is available through
#' ConScape, but should be interpreted cautiously because eigenvalue summaries
#' are intrinsically global.
#'
#' @param wrt Character vector describing the local perturbation target.
#'   Supported values are `"Q"` (habitat quality), `"A"` (affinity),
#'   `"C"` (cost), `"A&C=f(A)"` (linked affinity-to-cost change), and
#'   `"C&A=f(C)"` (linked cost-to-affinity change).
#' @param landscape_measure Summary of the landscape matrix: `"sum"` or
#'   `"eigenanalysis"`.
#' @param unitless Logical. If `TRUE`, return proportional sensitivities
#'   (elasticities). If `FALSE`, return raw unit-change sensitivities.
#' @param method `"analytical"` uses `ConScape.sensitivity`; `"simulation"`
#'   uses `ConScape.sensitivity_simulation`.
#' @param one_out_of Integer subsampling interval used only for simulation.
#' @param diagvalue Optional diagonal value passed to ConScape. Use `NULL` for
#'   ConScape's default.
#' @param target_equal_source Logical passed to ConScape. The current
#'   sensitivity implementation expects target quality to equal source quality.
#' @param require_landmark_one Logical. If `TRUE` (default),
#'   [run_conscape()] requires `landmark = 1L` for sensitivity runs so that
#'   ConScape's target qualities are not changed by coarse graining.
#'
#' @return An object of class `"ConScapeSensitivitySpec"`.
#' @export
conscape_sensitivity <- function(wrt = c("Q", "A&C=f(A)"),
                                 landscape_measure = c("sum", "eigenanalysis"),
                                 unitless = TRUE,
                                 method = c("analytical", "simulation"),
                                 one_out_of = 1L,
                                 diagvalue = NULL,
                                 target_equal_source = TRUE,
                                 require_landmark_one = TRUE) {
  allowed_wrt <- c("A", "C", "Q", "C&A=f(C)", "A&C=f(A)")
  if (missing(wrt)) {
    wrt <- c("Q", "A&C=f(A)")
  }
  if (!is.character(wrt) || length(wrt) == 0L) {
    stop("wrt must be a non-empty character vector.", call. = FALSE)
  }
  bad_wrt <- setdiff(wrt, allowed_wrt)
  if (length(bad_wrt) > 0L) {
    stop(
      "Unsupported sensitivity wrt value(s): ",
      paste(bad_wrt, collapse = ", "),
      call. = FALSE
    )
  }

  landscape_measure <- match.arg(landscape_measure)
  method <- match.arg(method)

  if (!is.logical(unitless) || length(unitless) != 1L || is.na(unitless)) {
    stop("unitless must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(one_out_of) || length(one_out_of) != 1L ||
      is.na(one_out_of) || one_out_of < 1) {
    stop("one_out_of must be a positive integer.", call. = FALSE)
  }
  if (!is.null(diagvalue) &&
      (!is.numeric(diagvalue) || length(diagvalue) != 1L || is.na(diagvalue))) {
    stop("diagvalue must be NULL or a single numeric value.", call. = FALSE)
  }

  out <- list(
    wrt = unique(wrt),
    landscape_measure = landscape_measure,
    unitless = unitless,
    method = method,
    one_out_of = as.integer(one_out_of),
    diagvalue = diagvalue,
    target_equal_source = isTRUE(target_equal_source),
    require_landmark_one = isTRUE(require_landmark_one)
  )
  class(out) <- "ConScapeSensitivitySpec"
  out
}

normalize_conscape_sensitivity <- function(x) {
  if (is.null(x) || identical(x, FALSE)) {
    return(NULL)
  }
  if (identical(x, TRUE)) {
    return(conscape_sensitivity())
  }
  if (inherits(x, "ConScapeSensitivitySpec")) {
    return(x)
  }
  if (is.character(x)) {
    return(conscape_sensitivity(wrt = x))
  }
  stop(
    "sensitivity must be NULL, FALSE, TRUE, a character vector, or a ConScapeSensitivitySpec.",
    call. = FALSE
  )
}

normalize_conscape_metrics <- function(metrics) {
  if (is.null(metrics)) {
    return(character())
  }
  if (!is.character(metrics)) {
    stop("metrics must be a character vector or NULL.", call. = FALSE)
  }

  aliases <- c(
    fcon = "connected_habitat",
    functional_connectivity = "connected_habitat",
    connected_habitat = "connected_habitat",
    btwn = "betweenness_kweighted",
    betweenness = "betweenness_kweighted",
    betweenness_kweighted = "betweenness_kweighted",
    betweenness_qweighted = "betweenness_qweighted",
    criticality = "criticality"
  )

  normalized <- unname(aliases[metrics])
  missing_alias <- is.na(normalized)
  if (any(missing_alias)) {
    stop(
      "Unsupported metric(s): ",
      paste(metrics[missing_alias], collapse = ", "),
      call. = FALSE
    )
  }
  unique(normalized)
}

normalize_conscape_connectivity_function <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop("connectivity_function must be a single character value.", call. = FALSE)
  }
  aliases <- c(
    expected_cost = "expected_cost",
    expected = "expected_cost",
    least_cost_distance = "least_cost_distance",
    least_cost = "least_cost_distance",
    free_energy_distance = "free_energy_distance",
    free_energy = "free_energy_distance",
    survival_probability = "survival_probability",
    survival = "survival_probability",
    power_mean_proximity = "power_mean_proximity",
    power_mean = "power_mean_proximity"
  )
  key <- tolower(trimws(x))
  out <- unname(aliases[key])
  if (is.na(out)) {
    stop("Unsupported connectivity_function: ", x, call. = FALSE)
  }
  out
}

normalize_conscape_cost_function <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop("cost_function must be a single character value.", call. = FALSE)
  }
  key <- tolower(trimws(x))
  aliases <- c(
    minuslog = "minuslog",
    "-log" = "minuslog",
    "x -> -log(x)" = "minuslog",
    inverse = "inverse",
    inv = "inverse",
    "x -> inv(x)" = "inverse",
    odds_against = "odds_against",
    odds_for = "odds_for",
    expminus = "expminus"
  )
  out <- unname(aliases[key])
  if (!is.na(out)) {
    return(out)
  }
  if (grepl("->", x, fixed = TRUE)) {
    return(x)
  }
  stop("Unsupported cost_function: ", x, call. = FALSE)
}

conscape_cost_function_julia <- function(cost_function, adjacency_expr) {
  cost_function <- normalize_conscape_cost_function(cost_function)
  switch(
    cost_function,
    minuslog = "ConScape.MinusLog()",
    inverse = "ConScape.Inv()",
    odds_against = "ConScape.OddsAgainst()",
    odds_for = "ConScape.OddsFor()",
    expminus = "ConScape.ExpMinus()",
    paste0("ConScape.mapnz(", cost_function, ", ", adjacency_expr, ")")
  )
}

conscape_metric_specs <- function() {
  list(
    connected_habitat = list(
      id = "connected_habitat",
      dir = "fcon",
      prefix = "fcon",
      layer = "fcon"
    ),
    betweenness_kweighted = list(
      id = "betweenness_kweighted",
      dir = "btwn",
      prefix = "betweenness",
      layer = "btwn"
    ),
    betweenness_qweighted = list(
      id = "betweenness_qweighted",
      dir = "btwn_qweighted",
      prefix = "betweenness_qweighted",
      layer = "btwn_qweighted"
    ),
    criticality = list(
      id = "criticality",
      dir = "criticality",
      prefix = "criticality",
      layer = "criticality"
    )
  )
}

conscape_sensitivity_id <- function(wrt, unitless) {
  base <- switch(
    wrt,
    Q = "quality",
    A = "affinity",
    C = "cost",
    "A&C=f(A)" = "affinity_cost_linked",
    "C&A=f(C)" = "cost_affinity_linked"
  )
  paste(if (unitless) "elasticity" else "sensitivity", base, sep = "_")
}

conscape_output_specs <- function(metrics, sensitivity = NULL) {
  metrics <- normalize_conscape_metrics(metrics)
  metric_catalog <- conscape_metric_specs()
  specs <- metric_catalog[metrics]

  sensitivity <- normalize_conscape_sensitivity(sensitivity)
  if (!is.null(sensitivity)) {
    sens_specs <- lapply(sensitivity$wrt, function(wrt) {
      id <- conscape_sensitivity_id(wrt, sensitivity$unitless)
      list(id = id, dir = id, prefix = id, layer = id, wrt = wrt)
    })
    names(sens_specs) <- vapply(sens_specs, `[[`, character(1), "id")
    specs <- c(specs, sens_specs)
  }

  if (length(specs) == 0L) {
    stop("At least one metric or sensitivity output must be requested.", call. = FALSE)
  }
  specs
}
