#' Assess ConScape prep workload
#'
#' @description
#' Summarize tile-level workload proxies from one or more objects created by
#' [conscape_prep()]. This is useful for comparing tile designs before
#' launching a costly Julia run.
#'
#' @param ... One or more `"ConScapeRtools_prep"` objects. Named arguments are
#'   used as scenario labels.
#' @param reference Integer or character identifying the reference scenario for
#'   reduction columns. Defaults to the first supplied prep object.
#'
#' @details
#' ConScape randomized shortest path calculations scale strongly with the
#' number of usable source graph cells and target quality cells in each tile.
#' The dominant cost per tile is solving `Z = (I - W) \ I_target` where the
#' number of columns of `Z` equals the number of positive target cells. Two
#' workload proxies are reported:
#'
#' * `estimated_source_target_cells` summarizes the work the ConScape backend
#'   actually performs per prep object. For `target_mode = "center"` (the new
#'   default), each tile's backend target raster is restricted to its center
#'   cells, so this column uses `n_positive_move * n_positive_output_target`.
#'   For `target_mode = "full"` (the legacy mode), the full buffered target
#'   raster is solved, so this column uses `n_positive_move * n_positive_target`.
#' * `estimated_output_source_target_cells` always uses
#'   `n_positive_move * n_positive_output_target`. It reflects the work that
#'   *would* be required if center-target windowing were used, regardless of
#'   the actual target mode of the prep.
#'
#' Likewise, `positive_target_cells` reports the count of positive cells in the
#' buffered tile's target raster (pre-restriction), and
#' `positive_output_target_cells` reports the count of positive target cells
#' that are written to disk and consumed by the backend. For
#' `target_mode = "center"`, the two differ; for `target_mode = "full"`, they
#' are equal.
#'
#' @return
#' A data frame with one row per prep object. Reduction columns are proportions
#' relative to `reference`, so `0.25` means a 25 percent reduction.
#'
#' @export
#' @examples
#' \dontrun{
#' classic <- conscape_prep(
#'   tile_d = 5000,
#'   tile_trim = 1000,
#'   r_target = habitat,
#'   r_mov = affinity,
#'   r_src = habitat,
#'   clear_dir = TRUE
#' )
#'
#' centered <- conscape_prep(
#'   centersize = 100,
#'   buffer = 20,
#'   r_target = habitat,
#'   r_mov = affinity,
#'   r_src = habitat,
#'   clear_dir = TRUE
#' )
#'
#' conscape_efficiency_assessment(classic = classic, centered = centered)
#' }
conscape_efficiency_assessment <- function(..., reference = 1L) {
  preps <- list(...)
  if (length(preps) == 0L) {
    stop("Supply at least one ConScapeRtools_prep object.", call. = FALSE)
  }

  labels <- names(preps)
  missing_labels <- !nzchar(labels)
  labels[missing_labels] <- paste0("prep_", which(missing_labels))

  summaries <- Map(function(prep, label) {
    if (!inherits(prep, "ConScapeRtools_prep")) {
      stop("All inputs must be ConScapeRtools_prep objects.", call. = FALSE)
    }
    diagnostics <- prep$tile_diagnostics
    if (is.null(diagnostics)) {
      diagnostics <- prep$tile_validation$diagnostics
    }
    if (is.null(diagnostics) || !is.data.frame(diagnostics)) {
      stop(
        "Prep object does not include tile diagnostics. Recreate it with the current conscape_prep().",
        call. = FALSE
      )
    }

    required <- c(
      "ncell",
      "n_positive_target",
      "n_positive_source",
      "n_positive_move",
      "n_finite_target",
      "n_finite_source",
      "n_finite_move"
    )
    missing <- setdiff(required, names(diagnostics))
    if (length(missing) > 0L) {
      stop("Prep diagnostics are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    if (!("n_positive_output_target" %in% names(diagnostics))) {
      diagnostics$n_positive_output_target <- diagnostics$n_positive_target
    }
    if (!("n_finite_output_target" %in% names(diagnostics))) {
      diagnostics$n_finite_output_target <- diagnostics$n_finite_target
    }

    # Choose which target count drives the backend-work proxy: with the
    # corrected center mode the backend solves Z over center-only targets, so
    # the proxy should use n_positive_output_target there.
    target_mode_for_work <- if (is.null(prep$target_mode)) "full" else prep$target_mode
    backend_target_count <- if (identical(target_mode_for_work, "center")) {
      diagnostics$n_positive_output_target
    } else {
      diagnostics$n_positive_target
    }
    estimated_source_target_cells <- sum(
      diagnostics$n_positive_move * backend_target_count,
      na.rm = TRUE
    )
    estimated_output_source_target_cells <- sum(
      diagnostics$n_positive_move * diagnostics$n_positive_output_target,
      na.rm = TRUE
    )
    data.frame(
      scenario = label,
      target_mode = if (is.null(prep$target_mode)) NA_character_ else prep$target_mode,
      n_tiles = length(prep$tile_num),
      landmark = prep$landmark,
      tile_trim = prep$tile_trim,
      centersize = if (is.null(prep$centersize)) NA_integer_ else prep$centersize,
      buffer = if (is.null(prep$buffer)) NA_integer_ else prep$buffer,
      tile_cells_total = sum(diagnostics$ncell, na.rm = TRUE),
      positive_target_cells = sum(diagnostics$n_positive_target, na.rm = TRUE),
      positive_output_target_cells = sum(diagnostics$n_positive_output_target, na.rm = TRUE),
      positive_source_cells = sum(diagnostics$n_positive_source, na.rm = TRUE),
      positive_movement_cells = sum(diagnostics$n_positive_move, na.rm = TRUE),
      finite_target_cells = sum(diagnostics$n_finite_target, na.rm = TRUE),
      finite_output_target_cells = sum(diagnostics$n_finite_output_target, na.rm = TRUE),
      finite_source_cells = sum(diagnostics$n_finite_source, na.rm = TRUE),
      finite_movement_cells = sum(diagnostics$n_finite_move, na.rm = TRUE),
      max_tile_positive_target_cells = max(diagnostics$n_positive_target, na.rm = TRUE),
      max_tile_positive_output_target_cells = max(diagnostics$n_positive_output_target, na.rm = TRUE),
      max_tile_positive_movement_cells = max(diagnostics$n_positive_move, na.rm = TRUE),
      estimated_source_target_cells = estimated_source_target_cells,
      estimated_output_source_target_cells = estimated_output_source_target_cells,
      stringsAsFactors = FALSE
    )
  }, preps, labels)

  out <- do.call(rbind, summaries)
  rownames(out) <- NULL

  reference_idx <- if (is.character(reference)) {
    match(reference, out$scenario)
  } else {
    as.integer(reference)
  }
  if (length(reference_idx) != 1L || is.na(reference_idx) ||
      reference_idx < 1L || reference_idx > nrow(out)) {
    stop("reference must identify one supplied prep object.", call. = FALSE)
  }

  ref_targets <- out$positive_target_cells[reference_idx]
  ref_output_targets <- out$positive_output_target_cells[reference_idx]
  ref_work <- out$estimated_source_target_cells[reference_idx]
  ref_output_work <- out$estimated_output_source_target_cells[reference_idx]
  out$target_cell_reduction <- if (ref_targets > 0) {
    1 - out$positive_target_cells / ref_targets
  } else {
    NA_real_
  }
  out$output_target_cell_reduction <- if (ref_output_targets > 0) {
    1 - out$positive_output_target_cells / ref_output_targets
  } else {
    NA_real_
  }
  out$source_target_work_reduction <- if (ref_work > 0) {
    1 - out$estimated_source_target_cells / ref_work
  } else {
    NA_real_
  }
  out$output_source_target_work_reduction <- if (ref_output_work > 0) {
    1 - out$estimated_output_source_target_cells / ref_output_work
  } else {
    NA_real_
  }

  class(out) <- c("ConScapeEfficiencyAssessment", class(out))
  out
}
