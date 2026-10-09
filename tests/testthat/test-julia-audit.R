test_that("Julia kernels and distributed scheduling satisfy regression contracts", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION"), "true"),
               "Set RUN_CONSCAPERTOOLS_INTEGRATION=true for Julia regression tests.")
  jl_home <- Sys.getenv("CONSCAPERTOOLS_JL_HOME", "")
  skip_if_not(nzchar(jl_home), "Set CONSCAPERTOOLS_JL_HOME.")
  exe <- ConScapeRtools:::conscape_julia_exe(jl_home)
  # test_path works from both a source checkout and R CMD check's test directory.
  package_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  if (!file.exists(file.path(package_root, "inst", "extdata", "conscape.jl"))) {
    package_root <- system.file(package = "ConScapeRtools")
    # Installed packages store extdata directly; these source-only driver tests
    # are also run explicitly in developer validation.
    skip("Julia audit drivers require a source checkout.")
  }
  for (script in c("audit_kernels.jl", "audit_scheduler.jl")) {
    output <- system2(exe, c("--startup-file=no", shQuote(testthat::test_path("..", "julia", script)),
                            shQuote(package_root)), stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    expect_true(is.null(status) || identical(status, 0L), info = paste(output, collapse = "\n"))
  }
})

test_that("real calibration meets its sampled threshold and releases Julia proxies", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION"), "true"),
               "Set RUN_CONSCAPERTOOLS_INTEGRATION=true for Julia regression tests.")
  jl_home <- Sys.getenv("CONSCAPERTOOLS_JL_HOME", "")
  skip_if_not(nzchar(jl_home), "Set CONSCAPERTOOLS_JL_HOME.")
  r <- make_test_raster(n = 5, vals = 0.5, crs = "EPSG:3857")
  design <- tile_design(r_mov = r, r_target = r, r_source = r, max_d = 2,
                        landmark = 1L, jl_home = jl_home)
  expect_lte(design$diagnostics$proximity_at_requested_trim, design$trim_threshold)
  expect_lte(design$diagnostics$proximity_at_effective_trim, design$trim_threshold)
  expect_false(conscape_julia_status())
})

test_that("all classic execution modes agree for asymmetric qualities and edge tiles", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RUN_CONSCAPERTOOLS_INTEGRATION"), "true"),
               "Set RUN_CONSCAPERTOOLS_INTEGRATION=true for Julia regression tests.")
  jl_home <- Sys.getenv("CONSCAPERTOOLS_JL_HOME", "")
  skip_if_not(nzchar(jl_home), "Set CONSCAPERTOOLS_JL_HOME.")
  mov <- terra::rast(nrows = 5, ncols = 7, xmin = 0, xmax = 7,
                     ymin = 0, ymax = 5, vals = seq(0.2, 0.8, length.out = 35),
                     crs = "EPSG:3857")
  source <- target <- mov
  terra::values(source) <- seq_len(35) / 35
  terra::values(target) <- rep(c(0, 0.3, 0.7, 0), length.out = 35)
  prep <- conscape_prep(centersize = 3, buffer = 7, landmark = 1L,
                        r_target = target, r_src = source, r_mov = mov, progress = FALSE)
  reference <- run_conscape(target_qualities = target, source_qualities = source,
                            affinities = mov, out_dir = tempfile("mode-reference-"),
                            landmark = 1L, jl_home = jl_home, distance_scale = 3,
                            progress = FALSE)
  modes <- list(serial = list(), threaded = list(parallel = TRUE, workers = 2L),
                distributed = list(parallel = TRUE, distributed = TRUE, workers = 2L),
                future = list(parallel = TRUE, parallel_R = TRUE, workers = 2L))
  for (mode in names(modes)) {
    out <- do.call(run_conscape, c(list(conscape_prep = prep, jl_home = jl_home,
                                        out_dir = tempfile(paste0(mode, "-")), distance_scale = 3,
                                        progress = FALSE), modes[[mode]]))
    expect_true(all(out$diagnostics$output_validation$ok), info = mode)
    for (metric in c("fcon", "btwn")) {
      expect_equal(terra::values(out[[metric]], mat = FALSE),
                   terra::values(reference[[metric]], mat = FALSE),
                   tolerance = 1e-6, info = paste(mode, metric))
    }
  }
  # A reused Julia session can contain caller-owned distributed workers.
  # This run must neither load our helpers there nor remove those workers.
  conscape_julia_start(jl_home, quiet = TRUE)
  on.exit(conscape_julia_stop(), add = TRUE)
  retained_worker <- JuliaConnectoR::juliaEval("using Distributed; only(addprocs(1))")
  out <- run_conscape(conscape_prep = prep, jl_home = jl_home,
                       out_dir = tempfile("existing-worker-"), distance_scale = 3,
                       progress = FALSE, parallel = TRUE, distributed = TRUE,
                       workers = 1L, stop_julia = FALSE)
  # Julia arrays carry JLDIM metadata; compare the worker IDs, not that attribute.
  expect_equal(as.integer(JuliaConnectoR::juliaEval("workers()")),
               as.integer(retained_worker))
  expect_false(JuliaConnectoR::juliaEval(paste0(
    "remotecall_fetch(Core.eval, ", retained_worker,
    ", Main, :(isdefined(Main, :conscape)))")))
  expect_false(JuliaConnectoR::juliaEval(paste0(
    "remotecall_fetch(Core.eval, ", retained_worker,
    ", Main, :(isdefined(Main, :safe_conscape)))")))
  expect_true(all(out$diagnostics$output_validation$ok))
})
