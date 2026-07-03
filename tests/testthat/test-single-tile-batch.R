## Regression coverage for a single-tile MethodError in the Julia-native
## parallel strategies (backend = "conscape_dev" is unaffected; this is the
## STABLE backend's julia_threads / julia_distributed batch path).
##
## Root cause: conscape_batch()/conscape_batch_distributed() in
## inst/extdata/conscape_batch.jl and conscape_batch_distributed.jl
## originally declared `r_targets::Vector{String}` (and r_sources, r_res).
## R does not distinguish a length-1 character vector from a scalar, so
## JuliaConnectoR auto-unboxes a single-tile design's one-element filename
## list into a bare Julia `String` instead of a `Vector{String}`. The strict
## type annotation then rejected the call with a MethodError, which fired
## 100% of the time whenever a tile design collapsed to exactly one tile
## (e.g., a landscape smaller than the calibrated minimum tile width) and
## the caller used `parallel = TRUE, parallel_R = FALSE` (julia_threads or
## julia_distributed). serial and future_R (R-side parallelism) never hit
## this because they call the per-tile `conscape()` Julia function once per
## tile with scalar filename arguments, which matches its scalar signature.
##
## Fix: relax the signatures to accept either form and normalize with
## `_as_string_vector()` (already defined in conscape.jl and used elsewhere
## in this codebase for the identical R-to-Julia scalar/vector ambiguity).

test_that("conscape_batch.jl no longer requires a strict Vector{String}", {
  # Static regression guard, no Julia required. Confirms the exact buggy
  # token is gone and the fix's normalization call is present, for both
  # the threaded and distributed batch scripts.
  batch_file <- system.file("extdata", "conscape_batch.jl", package = "ConScapeRtools")
  distributed_file <- system.file("extdata", "conscape_batch_distributed.jl", package = "ConScapeRtools")
  batch <- paste(readLines(batch_file, warn = FALSE), collapse = "\n")
  distributed <- paste(readLines(distributed_file, warn = FALSE), collapse = "\n")

  for (txt in list(batch = batch, distributed = distributed)) {
    expect_false(grepl("r_targets::Vector\\{String\\}", txt, fixed = FALSE))
    expect_false(grepl("r_sources::Vector\\{String\\}", txt, fixed = FALSE))
    expect_false(grepl("r_res::Vector\\{String\\}", txt, fixed = FALSE))
    expect_match(txt, "r_targets = _as_string_vector\\(r_targets\\)")
    expect_match(txt, "r_sources = _as_string_vector\\(r_sources\\)")
    expect_match(txt, "r_res\\s*= _as_string_vector\\(r_res\\)")
  }
})

integration_jl_home <- function() {
  candidates <- c(
    Sys.getenv("CONSCAPERTOOLS_JL_HOME", ""),
    Sys.getenv("JULIA_BINDIR", "")
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

single_tile_prep <- function(asc_dir) {
  # A 20 x 20 raster with tile_d spanning the whole extent collapses to
  # exactly one tile, which is the precondition that triggered the bug.
  r <- make_test_raster(n = 20, vals = 1)
  prep <- conscape_prep(
    tile_d    = 20,
    tile_trim = 2,
    asc_dir   = asc_dir,
    r_target  = r,
    r_src     = r,
    r_mov     = r,
    landmark  = 5L,
    clear_dir = TRUE,
    progress  = FALSE
  )
  expect_length(prep$tile_num, 1L)
  prep
}

test_that("julia_threads succeeds on a single-tile design (integration)", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION", ""), "true"),
    "Set RUN_CONSCAPERTOOLS_INTEGRATION=true to run Julia-backed ConScape checks."
  )
  jl_home <- integration_jl_home()
  skip_if_not(nzchar(jl_home), "Set CONSCAPERTOOLS_JL_HOME or JULIA_BINDIR.")

  prep <- single_tile_prep(file.path(tempdir(), "single-tile-threads-prep"))

  expect_no_error(
    out <- run_conscape(
      conscape_prep  = prep,
      out_dir        = file.path(tempdir(), "single-tile-threads-run"),
      theta          = 0.1,
      distance_scale = 10,
      jl_home        = jl_home,
      parallel       = TRUE,
      parallel_R     = FALSE,
      distributed    = FALSE,
      workers        = 2,
      progress       = FALSE
    )
  )
  expect_s3_class(out, "ConScapeResults")
})

test_that("julia_distributed succeeds on a single-tile design (integration)", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION", ""), "true"),
    "Set RUN_CONSCAPERTOOLS_INTEGRATION=true to run Julia-backed ConScape checks."
  )
  jl_home <- integration_jl_home()
  skip_if_not(nzchar(jl_home), "Set CONSCAPERTOOLS_JL_HOME or JULIA_BINDIR.")

  prep <- single_tile_prep(file.path(tempdir(), "single-tile-distributed-prep"))

  expect_no_error(
    out <- run_conscape(
      conscape_prep  = prep,
      out_dir        = file.path(tempdir(), "single-tile-distributed-run"),
      theta          = 0.1,
      distance_scale = 10,
      jl_home        = jl_home,
      parallel       = TRUE,
      parallel_R     = FALSE,
      distributed    = TRUE,
      workers        = 1,
      progress       = FALSE
    )
  )
  expect_s3_class(out, "ConScapeResults")
})
