## Functions from `ConScapeR` package

clear_juliaconnector_finalized_refs <- function() {
  ns <- tryCatch(asNamespace("JuliaConnectoR"), error = function(e) NULL)
  if (is.null(ns) || !exists("pkgLocal", envir = ns, inherits = FALSE)) {
    return(invisible(FALSE))
  }

  pkg_local <- get("pkgLocal", envir = ns, inherits = FALSE)
  if (!is.environment(pkg_local)) {
    return(invisible(FALSE))
  }

  pkg_local$finalizedRefs <- NULL
  invisible(TRUE)
}

juliaCall_conscape <- function(...) {
  clear_juliaconnector_finalized_refs()
  tryCatch(
    juliaCall(...),
    error = function(e) {
      clear_juliaconnector_finalized_refs()
      stop(e)
    }
  )
}

stop_conscape_julia <- function() {
  invisible(gc())
  try(stopJulia(), silent = TRUE)
  clear_juliaconnector_finalized_refs()
  invisible(NULL)
}

conscape_julia_exe <- function(jl_home) {
  candidates <- c(
    file.path(jl_home, "julia.exe"),
    file.path(jl_home, "julia"),
    jl_home
  )
  candidates[file.exists(candidates)][1]
}

default_conscape_dev_project <- function() {
  file.path(tools::R_user_dir("ConScapeRtools", which = "cache"),
            "julia", "conscape_dev")
}

run_julia_command <- function(command, args = character(), stdout = TRUE, stderr = TRUE) {
  system2(command, args = args, stdout = stdout, stderr = stderr)
}

prepare_conscape_dev_project <- function(jl_home,
                                         dev_project,
                                         install_dev_conscape,
                                         dev_conscape_rev,
                                         dev_conscape_url) {
  if (isFALSE(install_dev_conscape) && is.null(dev_project)) {
    return(NULL)
  }

  project <- if (is.null(dev_project)) default_conscape_dev_project() else dev_project
  if (isTRUE(install_dev_conscape)) {
    project <- conscape_dev_backend_setup(
      jl_home = jl_home,
      project = project,
      rev = dev_conscape_rev,
      url = dev_conscape_url,
      quiet = TRUE
    )
  } else if (!dir.exists(project)) {
    stop(
      "dev_project does not exist. Run conscape_dev_backend_setup() first, ",
      "or set install_dev_conscape = TRUE.",
      call. = FALSE
    )
  } else {
    project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  }

  current_project <- Sys.getenv("JULIA_PROJECT", unset = "")
  current_project <- if (nzchar(current_project)) {
    normalizePath(current_project, winslash = "/", mustWork = FALSE)
  } else {
    ""
  }
  if (!identical(current_project, project)) {
    if (conscape_julia_status()) {
      stop_conscape_julia()
    }
    Sys.setenv(JULIA_PROJECT = project)
  }
  project
}

#' Manage the Julia session used by ConScapeRtools
#'
#' @description
#' Start, reuse, inspect, or stop the Julia session used by ConScapeRtools.
#' These helpers are useful when running several Julia-backed operations in the
#' same R session and avoiding repeated Julia startup time.
#'
#' @param jl_home Path to the `bin` directory where the Julia executable
#'   is installed.
#' @param install_libraries Logical. If `TRUE`, install or update the required
#'   Julia libraries before importing them. Default is `FALSE`.
#' @param restart Logical. If `TRUE`, stop any existing Julia session before
#'   starting a new one. Use this after changing Julia paths, thread counts, or
#'   distributed-worker settings. Default is `FALSE`.
#' @param quiet Logical. If `TRUE`, suppress routine setup messages. Default is
#'   `FALSE`.
#' @param project Directory for a dedicated Julia project environment.
#'   Defaults to a package cache directory. Used by
#'   `conscape_dev_backend_setup()` only.
#' @param rev Git branch, tag, or commit for the ConScape dependency.
#'   `conscape_dev_backend_setup()` defaults to `"alg_efficiency"`;
#'   `conscape_sensitivity_setup()` defaults to `"sensitivity"`.
#' @param url Git URL for ConScape.jl.
#' @param force Logical. If `TRUE`, remove any existing ConScape dependency
#'   before adding the requested revision. For `conscape_sensitivity_setup()`
#'   this replaces a registered-release ConScape with the branch build (use
#'   when an older ConScape is pinned and `Pkg.add` will not move it).
#'
#' @return
#' `conscape_julia_start()` and `conscape_julia_stop()` return `invisible(NULL)`.
#' `conscape_julia_status()` returns a single logical value indicating whether
#' JuliaConnectoR currently reports an active Julia session.
#' `conscape_dev_backend_setup()` returns the normalized Julia project path.
#' `conscape_sensitivity_setup()` returns, invisibly, the installed ConScape
#' version string (e.g. `"0.3.0"`).
#'
#' @details
#' Reusing a Julia session is usually faster for repeated calls, but Julia
#' process settings such as `JULIA_NUM_THREADS` are fixed when Julia starts.
#' Restart Julia after changing those settings or after an error that may have
#' left Julia state uncertain.
#'
#' `conscape_sensitivity_setup()` installs a ConScape build that provides
#' `ConScape.sensitivity` and `ConScape.sensitivity_simulation` (and the
#' `criticality` and `betweenness_qweighted` metrics) into the **default**
#' Julia environment, the same environment [conscape_julia_start()] and
#' [run_conscape()] use. These functions live on the ConScape `sensitivity`
#' branch and are not part of the default registered ConScape release, so a
#' plain `Pkg.add("ConScape")` (what the package installs by default) does not
#' provide them. Call this helper once per machine before requesting
#' sensitivity outputs or the `criticality` / `betweenness_qweighted` metrics.
#' It is an explicit, opt-in install: ConScapeRtools never installs or swaps
#' Julia packages as a side effect of a normal [run_conscape()] call.
#'
#' @name conscape-julia
#' @aliases conscape_julia_start conscape_julia_stop conscape_julia_status
#'   conscape_dev_backend_setup conscape_sensitivity_setup
NULL

#' @rdname conscape-julia
#' @export
conscape_julia_start <- function(jl_home,
                                 install_libraries = FALSE,
                                 restart = FALSE,
                                 quiet = FALSE) {
  if (!is.character(jl_home) || length(jl_home) != 1L || is.na(jl_home) ||
      !nzchar(jl_home)) {
    stop("jl_home must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.logical(install_libraries) || length(install_libraries) != 1L ||
      is.na(install_libraries)) {
    stop("install_libraries must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(restart) || length(restart) != 1L || is.na(restart)) {
    stop("restart must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(restart)) {
    stop_conscape_julia()
  }

  clear_juliaconnector_finalized_refs()
  Sys.setenv(JULIA_BINDIR = jl_home)
  if (!suppressMessages(juliaSetupOk())) {
    stop("Check that the path to the Julia binary directory is correct", call. = FALSE)
  }

  setup_call <- function() ConScapeR_setup(jl_home, install_libraries = install_libraries)
  if (isTRUE(quiet)) {
    invisible(suppressMessages(setup_call()))
  } else {
    invisible(setup_call())
  }
}

#' @rdname conscape-julia
#' @export
conscape_julia_stop <- function() {
  stop_conscape_julia()
}

#' @rdname conscape-julia
#' @export
conscape_julia_status <- function() {
  clear_juliaconnector_finalized_refs()
  ns <- tryCatch(asNamespace("JuliaConnectoR"), error = function(e) NULL)
  if (is.null(ns) || !exists("pkgLocal", envir = ns, inherits = FALSE)) {
    return(FALSE)
  }

  pkg_local <- get("pkgLocal", envir = ns, inherits = FALSE)
  if (!is.environment(pkg_local) ||
      !exists("con", envir = pkg_local, inherits = FALSE) ||
      !exists("communicator", envir = pkg_local, inherits = FALSE)) {
    return(FALSE)
  }

  con <- get("con", envir = pkg_local, inherits = FALSE)
  communicator <- get("communicator", envir = pkg_local, inherits = FALSE)
  isTRUE(!is.null(communicator) &&
           inherits(con, "connection") &&
           isOpen(con))
}

#' @rdname conscape-julia
#' @export
conscape_dev_backend_setup <- function(jl_home,
                                       project = NULL,
                                       rev = "alg_efficiency",
                                       url = "https://github.com/ConScape/ConScape.jl",
                                       force = FALSE,
                                       quiet = FALSE) {
  if (!is.character(jl_home) || length(jl_home) != 1L || is.na(jl_home) ||
      !nzchar(jl_home)) {
    stop("jl_home must be a single non-empty character string.", call. = FALSE)
  }
  if (is.null(project)) {
    project <- default_conscape_dev_project()
  }
  if (!is.character(project) || length(project) != 1L || is.na(project) ||
      !nzchar(project)) {
    stop("project must be NULL or a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(rev) || length(rev) != 1L || is.na(rev) || !nzchar(rev)) {
    stop("rev must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(url) || length(url) != 1L || is.na(url) || !nzchar(url)) {
    stop("url must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.logical(force) || length(force) != 1L || is.na(force)) {
    stop("force must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE.", call. = FALSE)
  }

  julia <- conscape_julia_exe(jl_home)
  if (is.na(julia) || !nzchar(julia)) {
    stop("Check that the path to the Julia binary directory is correct", call. = FALSE)
  }

  dir.create(project, recursive = TRUE, showWarnings = FALSE)
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)

  setup_script <- tempfile("conscape_dev_setup_", fileext = ".jl")
  writeLines(c(
    "using Pkg",
    "project = ARGS[1]",
    "url = ARGS[2]",
    "rev = ARGS[3]",
    "force = ARGS[4] == \"true\"",
    "Pkg.activate(project; shared=false)",
    "if force",
    "    try",
    "        Pkg.rm(\"ConScape\")",
    "    catch",
    "    end",
    "end",
    "Pkg.add(Pkg.PackageSpec(url=url, rev=rev))",
    "Pkg.add(\"ArchGDAL\")",
    "Pkg.add(\"Rasters\")",
    "Pkg.instantiate()",
    "Pkg.precompile()",
    "using ConScape",
    "required = [:Problem, :WindowedProblem, :BatchProblem, :solve, :assess]",
    "missing = [String(s) for s in required if !isdefined(ConScape, s)]",
    "isempty(missing) || error(\"Installed ConScape is missing dev API symbols: \" * join(missing, \", \"))",
    "using ArchGDAL",
    "using Rasters",
    "println(\"ConScape dev backend project: \" * project)"
  ), setup_script, useBytes = TRUE)
  on.exit(unlink(setup_script, force = TRUE), add = TRUE)

  out <- run_julia_command(
    julia,
    args = c("--startup-file=no", setup_script, project, url, rev,
             if (isTRUE(force)) "true" else "false"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "Failed to prepare the ConScape development Julia project:\n",
      paste(out, collapse = "\n"),
      call. = FALSE
    )
  }
  if (isFALSE(quiet)) {
    message(paste(out, collapse = "\n"))
  }
  project
}

#' @rdname conscape-julia
#' @export
conscape_sensitivity_setup <- function(jl_home,
                                       rev = "sensitivity",
                                       url = "https://github.com/ConScape/ConScape.jl",
                                       force = FALSE,
                                       quiet = FALSE) {
  if (!is.character(jl_home) || length(jl_home) != 1L || is.na(jl_home) ||
      !nzchar(jl_home)) {
    stop("jl_home must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(rev) || length(rev) != 1L || is.na(rev) || !nzchar(rev)) {
    stop("rev must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(url) || length(url) != 1L || is.na(url) || !nzchar(url)) {
    stop("url must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.logical(force) || length(force) != 1L || is.na(force)) {
    stop("force must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE.", call. = FALSE)
  }

  julia <- conscape_julia_exe(jl_home)
  if (is.na(julia) || !nzchar(julia)) {
    stop("Check that the path to the Julia binary directory is correct", call. = FALSE)
  }

  # Install into the DEFAULT Julia environment (no Pkg.activate) so the
  # stable run_conscape() path, which does `using ConScape`, picks it up.
  # Shelling out to a fresh Julia process avoids mutating any live
  # JuliaConnectoR session's package state mid-run.
  setup_script <- tempfile("conscape_sensitivity_setup_", fileext = ".jl")
  writeLines(c(
    "using Pkg",
    "url = ARGS[1]",
    "rev = ARGS[2]",
    "force = ARGS[3] == \"true\"",
    "if force",
    "    try",
    "        Pkg.rm(\"ConScape\")",
    "    catch",
    "    end",
    "end",
    "Pkg.add(Pkg.PackageSpec(url=url, rev=rev))",
    "Pkg.precompile()",
    "using ConScape",
    "required = [:sensitivity, :sensitivity_simulation]",
    "missing = [String(s) for s in required if !isdefined(ConScape, s)]",
    "isempty(missing) || error(\"Installed ConScape is missing sensitivity API symbols: \" * join(missing, \", \"))",
    "println(\"ConScape sensitivity version: \" * string(pkgversion(ConScape)))"
  ), setup_script, useBytes = TRUE)
  on.exit(unlink(setup_script, force = TRUE), add = TRUE)

  out <- run_julia_command(
    julia,
    args = c("--startup-file=no", setup_script, url, rev,
             if (isTRUE(force)) "true" else "false"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "Failed to install the ConScape sensitivity build:\n",
      paste(out, collapse = "\n"),
      call. = FALSE
    )
  }
  if (isFALSE(quiet)) {
    message(paste(out, collapse = "\n"))
  }

  version_line <- grep("^ConScape sensitivity version: ", out, value = TRUE)
  version <- if (length(version_line)) {
    sub("^ConScape sensitivity version: ", "", version_line[[1]])
  } else {
    NA_character_
  }
  invisible(version)
}

julia_warn_or_error_expr <- function(expr) {
  paste0(
    "begin\n",
    "  using Logging\n",
    "  Logging.with_logger(Logging.ConsoleLogger(stderr, Logging.Warn)) do\n",
    expr,
    "\n  end\n",
    "end"
  )
}

juliaLet_warn_or_error <- function(expr, ...) {
  juliaLet(julia_warn_or_error_expr(expr), ...)
}

#' Launch the Julia installation from R and load the `ConScape` library in Julia
#'
#' This function starts a Julia session from R and imports the `ConScape` library
#' to be used from R. It assumes that Julia is already installed in the system
#' and the path to its executables is given as the `julia_path` argument.
#' The first time the function is ran, it is best to set `install_libraries = TRUE`
#' to install the `ConScape` library in Julia.
#'
#' @param julia_path `[character]` \cr The directory for the Julia bin, e.g.
#' "C:/Programs/Julia-1.9.3/bin".
#' @param install_libraries `[logical]` \cr If `FALSE` (default), Julia will be
#' launched and the required libraries will be loaded without installing them.
#' If `TRUE`, the libraries will be (re-)installed in Julia.
#'
#' @return NULL.
#' @keywords internal
#' @noRd
#' @importFrom JuliaConnectoR juliaLet
#'
ConScapeR_setup <- function(julia_path, install_libraries = FALSE) {

  clear_juliaconnector_finalized_refs()
  Sys.setenv(JULIA_BINDIR = julia_path)

  # List of required packages
  required_pkgs <- c("ConScape", "SparseArrays", "Statistics")

  # Run the check and install function
  check_and_install_julia_pkgs(required_pkgs)

  if (install_libraries){
    Pkg <- juliaImport("Pkg")
    juliaEval("Pkg.add(\"ConScape\")")
    juliaEval("Pkg.add(\"SparseArrays\")")
    juliaEval("Pkg.add(\"Statistics\")")
  }
  SA <- juliaImport("SparseArrays")
  CS <- juliaImport("ConScape")
  invisible(juliaEval("using Logging"))
}

#' Wrapper for the Grid function of `ConScape`
#'
#' Creates a Grid from affinities, sources, targets and costs.
#' Affinities, sources and targets can be provided either as
#' `SpatRaster` from the `terra` library or as `matrix`.
#' The costs can be provided as `SpatRaster` or `matrix`,
#' but also as a string describing a transformation from the affinity matrix
#' (e.g. `"x -> -log(x)"`).
#'
#' @param affinities `[SpatRaster, matrix]` \cr The affinities, represent the likelihood
#' of movement. They can be provided as a `SpatRaster` representing the affinity,
#' permeability, or resistance of each pixel to movement (in which case it is
#' internally transformed into a matrix), or as a matrix with affinities between
#' pairs of sources-targets directly.
#' @param sources `[SpatRaster, matrix]` \cr Map or matrix presenting the habitat
#' suitability/quality of source cells.
#' @param targets `[SpatRaster, matrix]` \cr Map or matrix presenting the habitat
#' suitability/quality of target cells.
#' @param costs `[SpatRaster, matrix, character]` \cr Map or matrix presenting the
#' cost of movement through each cell. Alternatively, a string describing a transformation
#' from the affinity matrix (e.g. `"x -> -log(x)"`)
#'
#' @return A `ConScape.Grid` object within Julia.
#' @keywords internal
#' @noRd
Grid <- function(affinities, sources, targets, costs) {
  # affinities
  if (inherits(affinities, "SpatRaster")) {
    affinities <- terra::as.matrix(affinities, wide = T)
  }
  # sources
  if (inherits(sources, "SpatRaster")) {
    sources <- terra::as.matrix(sources, wide = T)
  }
  # targets
  if (inherits(targets, "SpatRaster")) {
    targets <- terra::as.matrix(targets, wide = T)
  }
  # costs
  if (inherits(costs, "SpatRaster")) {
    costs <- terra::as.matrix(costs, wide = T)
  }

  # check NaNs
  if (inherits(costs, "matrix")) {
    nans <- is.nan(sources) | is.nan(targets) | is.nan(affinities) | is.nan(costs)
    costs[nans] <- NaN
  } else {
    nans <- is.nan(sources) | is.nan(targets) | is.nan(affinities)
  }
  sources[nans] <- NaN
  targets[nans] <- NaN
  affinities[nans] <- NaN

  # create Grid
  if (inherits(costs, "character")){
    cost_expr <- conscape_cost_function_julia(
      costs,
      "ConScape.graph_matrix_from_raster(affinities)"
    )
    g <- juliaLet_warn_or_error(
      paste0("ConScape.Grid(size(affinities)...,
              affinities=ConScape.graph_matrix_from_raster(affinities),
              source_qualities=sources,
              target_qualities=SparseArrays.sparse(targets),
              costs=", cost_expr, ")"),
      affinities=affinities, sources=sources, targets=targets)
  } else {
    g <- juliaLet_warn_or_error("ConScape.Grid(size(affinities)...,
                  affinities=ConScape.graph_matrix_from_raster(affinities),
                  source_qualities=sources,
                  target_qualities=SparseArrays.sparse(targets),
                  costs=ConScape.graph_matrix_from_raster(costs))",
                                  affinities=affinities, sources=sources, targets=targets, costs=costs)
  }
  return(g)
}


#' Convert a matrix to raster
#'
#' @param mat `[matrix]` \cr Matrix to be converted.
#' @param rast `[SpatRaster]` \cr Template raster, usually one of those used in the
#' `ConScapeR::Grid` function.
#'
#' @return `[SpatRaster]`
#' @keywords internal
#' @noRd
#'
mat2rast <- function(mat, rast){
  rast2 <- terra::rast(mat, extent = terra::ext(rast), crs = terra::crs(rast))
  return(rast2)
}

#' Convert a node vector to matrix
#'
#' @param vec `[matrix]` \cr Matrix to be converted
#' @param g `[Grid]` \cr Template raster, usually one of those used in the
#' `ConScapeR::Grid` function.
#'
#' @return `[matrix]`
#' @keywords internal
#' @noRd
#'
vec2mat <- function(vec, g){
  mat <- juliaLet_warn_or_error("ConScape._vec_to_grid(g, vec)", g=g, vec=vec)
  return(mat)
}

#' Wrapper for the GridRSP function of `ConScape`
#'
#' Creates a GridRSP from a ConScape.Grid and the randomness parameter. Theta closer to zero gives a more random movement,
#' as theta increases the RSP distribution will converge to the optimal or least-cost paths.
#' Note that the computation becomes numerically unstable as theta becomes really small or large.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param g `[Grid]` \cr The output from the `ConScapeR::Grid` function.
#' @param theta `[numeric]` \cr The randomness parameter theta. Lower `theta` values
#' represent a more random walk, and higher values a more least-cost path behavior.
#'
#' @return A `ConScape.GridRSP` object within Julia.
#' @keywords internal
#' @noRd
#'
GridRSP <- function(g, theta) {
  h <- juliaLet_warn_or_error("ConScape.GridRSP(g, \u03b8=theta)", g=g, theta=theta)
  return(h)
}

#' Wrapper for the `betweenness_qweighted` function of `ConScape`
#'
#' Computes the quality-weighted betweenness for a ConScape.GridRSP object.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the `ConScapeR::Grid` function.
#'
#' @return `matrix`
#' @keywords internal
#' @noRd
#'
betweenness_qweighted <- function(h) {
  betw <- juliaLet_warn_or_error("ConScape.betweenness_qweighted(h)", h=h)
  return(betw)
}

#' Wrapper for the `betweenness_kweighted` function of `ConScape`
#'
#' Computes the quality- and proximity-weighted betweenness for a ConScape.GridRSP object.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the `ConScapeR::Grid` function.
#' @param alpha `[numeric]` \cr The distance scaling for the exponential distance
#' to proximity transformation.
#'
#' @return `matrix`
#' @keywords internal
#' @noRd
#'
betweenness_kweighted <- function(h, alpha) {
  betw <- juliaLet_warn_or_error("ConScape.betweenness_kweighted(h, distance_transformation=x -> exp(-x*alpha))", h=h, alpha=alpha)
  return(betw)
}

coarse_grid <- function(g, land_mark){
  coarse <- juliaLet_warn_or_error("ConScape.coarse_graining(g, land_mark)",
                      g=g, land_mark = land_mark)
  return(coarse)
}

#' Wrapper for the `connected_habitat` function of `ConScape`
#'
#' Computes the habitat functionality for a ConScape.GridRSP object.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the `ConScapeR::Grid` function.
#' @param alpha `[numeric]` \cr The distance scaling for the exponential distance to
#' proximity transformation.
#'
#' @return A `matrix` with loca (pixel) measures of habitat functionality.
#' @keywords internal
#' @noRd
#'
connected_habitat <- function(h, alpha) {
  func <- juliaLet_warn_or_error("ConScape.connected_habitat(h, distance_transformation=x -> exp(-x*alpha))", h=h, alpha=alpha)
  return(func)
}

#' Wrapper for the `expected_cost` function of `ConScape`
#'
#' Computes the RSP expected cost between all sources and non-zero targets.
#' See Van Moorter et al. (2023, Methods in Ecology and Evolution) for more details.
#'
#' @param h `[GridRSP]` \cr The output from the `ConScapeR::Grid` function.
#'
#' @keywords internal
#' @noRd
#'
expected_cost <- function(h) {
  dists = juliaLet_warn_or_error("ConScape.expected_cost(h)", h=h)
  return(dists)
}

sensitivity <- function(h,
                        wrt = "Q",
                        alpha = 1,
                        landscape_measure = "sum",
                        unitless = TRUE,
                        method = "analytical",
                        diagvalue = NULL,
                        target_equal_source = TRUE,
                        one_out_of = 1L) {
  fun <- if (identical(method, "simulation")) {
    "sensitivity_simulation"
  } else {
    "sensitivity"
  }
  extra <- if (identical(method, "simulation")) {
    ", one_out_of=one_out_of"
  } else {
    ""
  }
  juliaLet_warn_or_error(
    paste0(
      "ConScape.", fun,
      "(h, distance_transformation=ConScape.ExpMinus(), \u03b1=alpha, ",
      "wrt=wrt, landscape_measure=landscape_measure, unitless=unitless, ",
      "diagvalue=diagvalue, target_equal_source=target_equal_source", extra, ")"
    ),
    h = h,
    wrt = wrt,
    alpha = alpha,
    landscape_measure = landscape_measure,
    unitless = unitless,
    diagvalue = diagvalue,
    target_equal_source = target_equal_source,
    one_out_of = as.integer(one_out_of)
  )
}





# Function to check and install packages
check_and_install_julia_pkgs <- function(pkgs) {

  # Check each package
  for (pkg in pkgs) {
    # Improved check that properly handles the catch clause
    is_installed <- juliaEval(paste0(
      'using Pkg; ',
      'try ',
      '    using ', pkg, '; ',
      '    true ',
      'catch e ',
      '    false ',
      'end'
    ))

    if (!is_installed) {
      message("Installing Julia package: ", pkg)
      tryCatch({
        juliaEval(paste0('using Pkg; Pkg.add("', pkg, '")'))
        message("Successfully installed: ", pkg)
      }, error = function(e) {
        warning("Failed to install package ", pkg, ": ", e$message)
      })
    } else {
      message("Julia package ", pkg, " is already installed")
    }
  }

  # Final verification that packages can be loaded
  message("\nVerifying package loading:")
  for (pkg in pkgs) {
    tryCatch({
      juliaEval(paste0('using ', pkg))
      message("Successfully loaded: ", pkg)
    }, error = function(e) {
      warning("Failed to load package ", pkg, ": ", e$message)
    })
  }
}
