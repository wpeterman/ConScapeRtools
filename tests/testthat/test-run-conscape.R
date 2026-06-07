mock_julia <- function() {
  testthat::local_mocked_bindings(
    stopJulia = function() invisible(TRUE),
    juliaSetupOk = function() TRUE,
    juliaEval = function(...) invisible(TRUE),
    ConScapeR_setup = function(...) invisible(TRUE),
    .env = asNamespace("ConScapeRtools")
  )
}

write_fake_conscape_outputs <- function(out_dir, source_file, iter = "") {
  r <- terra::rast(source_file)
  btwn <- r
  fcon <- r
  terra::values(btwn) <- seq_len(terra::ncell(btwn))
  terra::values(fcon) <- seq_len(terra::ncell(fcon)) * 2
  dir.create(file.path(out_dir, "btwn"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "fcon"), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(btwn, file.path(out_dir, "btwn", paste0("betweenness", iter, ".asc")),
                     overwrite = TRUE, NAflag = -9999)
  terra::writeRaster(fcon, file.path(out_dir, "fcon", paste0("fcon", iter, ".asc")),
                     overwrite = TRUE, NAflag = -9999)
  "done"
}

test_that("run_conscape rejects non-empty output directory when requested", {
  mock_julia()
  out_dir <- file.path(tempdir(), "run-non-empty")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("existing", file.path(out_dir, "existing.txt"))

  expect_error(
    run_conscape(
      out_dir = out_dir,
      hab_target = "missing",
      hab_src = "missing",
      mov_prob = "missing",
      clear_dir = FALSE,
      jl_home = "C:/Julia/bin"
    ),
    "Files currently exist"
  )
})

test_that("run_conscape recreates output directory immediately after clearing it", {
  mock_julia()
  testthat::local_mocked_bindings(
    juliaCall = function(...) stop("simulated Julia failure"),
    .env = asNamespace("ConScapeRtools")
  )
  out_dir <- file.path(tempdir(), "run-clear-recreate")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("existing", file.path(out_dir, "existing.txt"))
  r <- make_test_raster(n = 3, vals = 1)

  expect_error(
    run_conscape(
      out_dir = out_dir,
      hab_target = r,
      hab_src = r,
      mov_prob = r,
      jl_home = "C:/Julia/bin"
    ),
    "simulated Julia failure"
  )
  expect_true(dir.exists(out_dir))
  expect_false(file.exists(file.path(out_dir, "existing.txt")))
})

test_that("run_conscape validates single-raster geometry before running Julia", {
  mock_julia()
  target <- make_test_raster()
  src <- make_test_raster(n = 10)
  mov <- make_test_raster()

  expect_error(
    run_conscape(
      out_dir = file.path(tempdir(), "run-bad-geometry"),
      hab_target = target,
      hab_src = src,
      mov_prob = mov,
      jl_home = "C:/Julia/bin"
    ),
    "share extent"
  )
})

test_that("run_conscape can leave Julia open when requested", {
  stop_calls <- 0L
  testthat::local_mocked_bindings(
    stopJulia = function() {
      stop_calls <<- stop_calls + 1L
      invisible(TRUE)
    },
    juliaSetupOk = function() TRUE,
    juliaEval = function(...) invisible(TRUE),
    ConScapeR_setup = function(...) invisible(TRUE),
    .env = asNamespace("ConScapeRtools")
  )

  target <- make_test_raster()
  src <- make_test_raster(n = 10)
  mov <- make_test_raster()

  expect_error(
    run_conscape(
      out_dir = file.path(tempdir(), "run-bad-geometry-keep-julia"),
      hab_target = target,
      hab_src = src,
      mov_prob = mov,
      jl_home = "C:/Julia/bin",
      stop_julia = FALSE
    ),
    "share extent"
  )
  expect_equal(stop_calls, 0L)
})

test_that("run_conscape warns when reusing Julia with threaded settings", {
  mock_julia()
  testthat::local_mocked_bindings(
    juliaEval = function(expr, ...) {
      if (identical(expr, "Threads.nthreads()")) return(1L)
      invisible(TRUE)
    },
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir, hab_target, ...) {
      for (i in seq_along(hab_target)) {
        write_fake_conscape_outputs(
          out_dir,
          file.path(target_dir, hab_target[[i]]),
          iter = paste0("-iter_", i)
        )
      }
      invisible(TRUE)
    },
    conscape_julia_status = function() TRUE,
    .env = asNamespace("ConScapeRtools")
  )

  r <- make_test_raster(n = 10, vals = 1)
  prep <- conscape_prep(
    tile_d = 5,
    tile_trim = 2,
    asc_dir = file.path(tempdir(), "run-thread-reuse-prep"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )
  expect_warning(
    run_conscape(
      out_dir = file.path(tempdir(), "run-thread-reuse-warning"),
      conscape_prep = prep,
      jl_home = "C:/Julia/bin",
      parallel = TRUE,
      workers = 2,
      stop_julia = FALSE
    ),
    "may not honor the requested workers"
  )
})

test_that("run_conscape can run a mocked single-raster workflow and crop extension", {
  mock_julia()
  calls <- list()
  testthat::local_mocked_bindings(
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir,
                         r_target, r_source, r_res, landmark, theta, exp_d, NA_val, iter, ...) {
      calls[[length(calls) + 1L]] <<- list(
        name = name, r_target = r_target, r_source = r_source, r_res = r_res,
        landmark = landmark, theta = theta, exp_d = exp_d, iter = iter
      )
      write_fake_conscape_outputs(out_dir, file.path(target_dir, r_target), iter)
    },
    .env = asNamespace("ConScapeRtools")
  )

  r <- make_test_raster(n = 6, vals = 1, crs = "EPSG:4326")
  out <- run_conscape(
    out_dir = file.path(tempdir(), "run-single"),
    hab_target = r,
    hab_src = r,
    mov_prob = r,
    landmark = 3L,
    theta = 0.2,
    exp_d = 42,
    tile_trim = 2,
    jl_home = "C:/Julia/bin"
  )

  expect_s4_class(out, "SpatRaster")
  expect_equal(names(out), c("btwn", "fcon"))
  expect_equal(dim(out), c(6, 6, 2))
  expect_equal(terra::crs(out), terra::crs(r))
  expect_equal(calls[[1]]$name, "conscape")
  expect_equal(calls[[1]]$landmark, 3L)
  expect_equal(calls[[1]]$iter, "")
})

test_that("run_conscape can run mocked tiled serial workflow and mosaic results", {
  mock_julia()
  r <- make_test_raster(n = 10, vals = 1)
  prep <- conscape_prep(
    tile_d = 5,
    tile_trim = 2,
    asc_dir = file.path(tempdir(), "run-prep"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  testthat::local_mocked_bindings(
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir,
                         r_target, r_source, r_res, landmark, theta, exp_d, NA_val, iter, ...) {
      write_fake_conscape_outputs(out_dir, file.path(target_dir, r_target), iter)
    },
    .env = asNamespace("ConScapeRtools")
  )

  out <- run_conscape(
    conscape_prep = prep,
    out_dir = file.path(tempdir(), "run-tiled"),
    theta = 0.1,
    exp_d = 50,
    jl_home = "C:/Julia/bin",
    progress = FALSE
  )

  expect_s3_class(out, "ConScapeResults")
  expect_s4_class(out$btwn, "SpatRaster")
  expect_s4_class(out$fcon, "SpatRaster")
  expect_equal(dim(out$btwn), dim(r))
  expect_true(dir.exists(out$outdir_btwn))
  expect_true(dir.exists(out$outdir_fcon))
})

test_that("run_conscape handles a prepared single-tile workflow", {
  mock_julia()
  r <- make_test_raster(n = 10, vals = 1)
  prep <- conscape_prep(
    tile_d = 10,
    tile_trim = 2,
    asc_dir = file.path(tempdir(), "run-single-tile-prep"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  testthat::local_mocked_bindings(
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir,
                         r_target, r_source, r_res, landmark, theta, exp_d, NA_val, iter, ...) {
      write_fake_conscape_outputs(out_dir, file.path(target_dir, r_target), iter)
    },
    .env = asNamespace("ConScapeRtools")
  )

  out <- run_conscape(
    conscape_prep = prep,
    out_dir = file.path(tempdir(), "run-single-tile"),
    theta = 0.1,
    exp_d = 50,
    jl_home = "C:/Julia/bin",
    progress = FALSE
  )

  expect_s3_class(out, "ConScapeResults")
  expect_equal(dim(out$btwn), dim(r))
})

test_that("run_conscape warns when mocked threaded execution misses outputs", {
  mock_julia()
  r <- make_test_raster(n = 10, vals = 1)
  prep <- conscape_prep(
    tile_d = 5,
    tile_trim = 2,
    asc_dir = file.path(tempdir(), "run-parallel-prep"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  testthat::local_mocked_bindings(
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir, ...) {
      dir.create(file.path(out_dir, "btwn"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(out_dir, "fcon"), recursive = TRUE, showWarnings = FALSE)
      "ok"
    },
    .env = asNamespace("ConScapeRtools")
  )

  expect_warning(
    out <- run_conscape(
      conscape_prep = prep,
      out_dir = file.path(tempdir(), "run-parallel"),
      theta = 0.1,
      exp_d = 50,
      jl_home = "C:/Julia/bin",
      parallel = TRUE,
      workers = 2,
      mosaic = FALSE
    ),
    "Parallel execution failed"
  )
  expect_s3_class(out, "ConScapeResults")
  expect_true("diagnostics" %in% names(out))
  expect_true("output_validation" %in% names(out$diagnostics))
})

test_that("run_conscape validates requested parallel outputs instead of defaults", {
  mock_julia()
  r <- make_test_raster(n = 10, vals = 1)
  prep <- conscape_prep(
    tile_d = 5,
    tile_trim = 2,
    asc_dir = file.path(tempdir(), "run-parallel-criticality-prep"),
    r_target = r,
    r_mov = r,
    r_src = r,
    clear_dir = TRUE,
    landmark = 5L,
    progress = FALSE
  )

  testthat::local_mocked_bindings(
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir,
                         hab_target, hab_src, mov_prob, ...) {
      dir.create(file.path(out_dir, "criticality"), recursive = TRUE, showWarnings = FALSE)
      for (i in seq_along(hab_target)) {
        tile <- terra::rast(file.path(target_dir, hab_target[[i]]))
        terra::writeRaster(
          tile,
          file.path(out_dir, "criticality", paste0("criticality-r_", i, ".asc")),
          overwrite = TRUE,
          NAflag = -9999
        )
      }
      "ok"
    },
    .env = asNamespace("ConScapeRtools")
  )

  expect_warning(
    out <- run_conscape(
      conscape_prep = prep,
      out_dir = file.path(tempdir(), "run-parallel-criticality"),
      theta = 0.1,
      exp_d = 50,
      jl_home = "C:/Julia/bin",
      parallel = TRUE,
      workers = 2,
      metrics = "criticality",
      mosaic = FALSE
    ),
    NA
  )
  expect_s3_class(out, "ConScapeResults")
  expect_true("criticality" %in% names(out$outdirs))
  expect_false("btwn" %in% names(out$outdirs))
  expect_true(all(out$diagnostics$output_validation$ok))
})

test_that("run_conscape accepts ConScape-native slot aliases", {
  mock_julia()
  calls <- list()
  testthat::local_mocked_bindings(
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir,
                         r_target, r_source, r_res, landmark, theta, exp_d, NA_val, iter, ...) {
      extras <- list(...)
      calls[[length(calls) + 1L]] <<- extras
      write_fake_conscape_outputs(out_dir, file.path(target_dir, r_target), iter)
    },
    .env = asNamespace("ConScapeRtools")
  )

  r <- make_test_raster(n = 6, vals = 1)
  out <- run_conscape(
    out_dir = file.path(tempdir(), "run-aliases"),
    target_qualities = r,
    source_qualities = r,
    affinities = r,
    distance_scale = 30,
    cost_function = "minuslog",
    connectivity_function = "expected",
    metrics = c("btwn", "fcon"),
    landmark = 3L,
    theta = 0.2,
    jl_home = "C:/Julia/bin"
  )

  expect_s4_class(out, "SpatRaster")
  expect_equal(names(out), c("btwn", "fcon"))
  expect_equal(calls[[1]][[1]], c("betweenness_kweighted", "connected_habitat"))
  expect_equal(calls[[1]][[2]], "expected_cost")
  expect_equal(calls[[1]][[3]], "minuslog")
})

test_that("run_conscape forwards experimental ConScape dev modes", {
  calls <- list()
  testthat::local_mocked_bindings(
    run_conscape_dev_backend = function(...) {
      calls[[length(calls) + 1L]] <<- list(...)
      out <- make_test_raster(n = 3, vals = 1)
      names(out) <- "btwn"
      out
    },
    .env = asNamespace("ConScapeRtools")
  )

  r <- make_test_raster(n = 6, vals = 1)
  out <- run_conscape(
    out_dir = file.path(tempdir(), "run-dev-batch"),
    target_qualities = r,
    source_qualities = r,
    affinities = r,
    backend = "conscape_dev",
    dev_mode = "batch",
    centersize = 3,
    buffer = 1,
    batch_grain = 2,
    batch_ext = "tif",
    dev_project = "C:/Julia/conscape-dev-project",
    install_dev_conscape = TRUE,
    dev_conscape_rev = "alg_efficiency",
    dev_conscape_url = "https://github.com/ConScape/ConScape.jl",
    landmark = 1L,
    jl_home = "C:/Julia/bin"
  )

  expect_s4_class(out, "SpatRaster")
  expect_equal(calls[[1]]$dev_mode, "batch")
  expect_equal(calls[[1]]$batch_grain, 2)
  expect_equal(calls[[1]]$batch_ext, "tif")
  expect_equal(calls[[1]]$dev_project, "C:/Julia/conscape-dev-project")
  expect_true(calls[[1]]$install_dev_conscape)
  expect_equal(calls[[1]]$dev_conscape_rev, "alg_efficiency")
  expect_equal(calls[[1]]$dev_conscape_url, "https://github.com/ConScape/ConScape.jl")
  expect_equal(calls[[1]]$window_shape, "square")
})

test_that("run_conscape validates experimental batch options early", {
  testthat::local_mocked_bindings(
    run_conscape_dev_backend = function(...) {
      list(...)
    },
    .env = asNamespace("ConScapeRtools")
  )

  r <- make_test_raster(n = 3, vals = 1)
  expect_error(
    run_conscape(
      out_dir = file.path(tempdir(), "bad-dev-batch"),
      target_qualities = r,
      source_qualities = r,
      affinities = r,
      backend = "conscape_dev",
      dev_mode = "batch",
      landmark = 1L,
      centersize = 3,
      buffer = 1,
      batch_grain = 0,
      jl_home = "C:/Julia/bin"
    ),
    "batch_grain"
  )
})

test_that("run_conscape guards sensitivity against coarse graining by default", {
  mock_julia()
  r <- make_test_raster(n = 6, vals = 1)

  expect_error(
    run_conscape(
      out_dir = file.path(tempdir(), "run-sensitivity-guard"),
      target_qualities = r,
      source_qualities = r,
      affinities = r,
      landmark = 3L,
      sensitivity = conscape_sensitivity(wrt = "Q"),
      jl_home = "C:/Julia/bin"
    ),
    "landmark = 1L"
  )
})

test_that("run_conscape reports prep landmark requirement for sensitivity", {
  # Sensitivity requires target_mode = "full" (Fix 1: center mode is rejected
  # because it violates ConScape.sensitivity's target_equal_source = TRUE
  # precondition). Use full-mode prep here so this test exercises the
  # landmark-precondition check rather than the center-mode rejection.
  mock_julia()
  r <- make_test_raster(n = 6, vals = 1)
  prep <- conscape_prep(
    tile_d = 6,
    tile_trim = 0,
    asc_dir = file.path(tempdir(), "run-sensitivity-prep-guard"),
    r_target = r,
    r_src = r,
    r_mov = r,
    landmark = 3L,
    target_mode = "full",
    clear_dir = TRUE,
    progress = FALSE
  )

  expect_error(
    run_conscape(
      conscape_prep = prep,
      out_dir = file.path(tempdir(), "run-sensitivity-prep-guard-out"),
      landmark = 1L,
      sensitivity = conscape_sensitivity(wrt = "Q"),
      jl_home = "C:/Julia/bin"
    ),
    "create that object with landmark = 1L"
  )
})

test_that("run_conscape rejects sensitivity with center-mode prep (Fix 1)", {
  # Center-mode prep + sensitivity request must error early with a clear
  # workaround message, before any Julia start or out_dir clearing.
  mock_julia()
  r <- make_test_raster(n = 6, vals = 1)
  prep <- conscape_prep(
    tile_d = 6,
    tile_trim = 0,
    asc_dir = file.path(tempdir(), "run-sensitivity-center-guard"),
    r_target = r,
    r_src = r,
    r_mov = r,
    landmark = 1L,
    target_mode = "center",
    clear_dir = TRUE,
    progress = FALSE
  )
  expect_identical(prep$target_mode, "center")

  expect_error(
    run_conscape(
      conscape_prep = prep,
      out_dir = file.path(tempdir(), "run-sensitivity-center-guard-out"),
      landmark = 1L,
      sensitivity = conscape_sensitivity(wrt = "Q"),
      jl_home = "C:/Julia/bin"
    ),
    "target_mode = \"center\""
  )
})
