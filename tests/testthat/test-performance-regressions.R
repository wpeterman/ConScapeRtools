test_that("directory inputs process every matched tile", {
  mock_julia()
  r <- make_test_raster(n = 4, vals = 0.5, crs = "EPSG:3857")
  prep <- conscape_prep(centersize = 2, buffer = 1, landmark = 1L,
                        r_target = r, r_mov = r, r_src = r, progress = FALSE)
  seen <- character()
  testthat::local_mocked_bindings(
    juliaCall = function(name, src_dir, mov_dir, target_dir, out_dir,
                         r_target, r_source, r_res, landmark, theta, exp_d, NA_val, iter, ...) {
      seen <<- c(seen, r_target)
      write_fake_conscape_outputs(out_dir, file.path(target_dir, r_target), iter)
    }, .package = "ConScapeRtools"
  )
  out <- run_conscape(target_qualities = prep$target, source_qualities = prep$src,
                      affinities = prep$mov, out_dir = tempfile("directory-run-"),
                      landmark = 1L, jl_home = "mock", progress = FALSE)
  expect_identical(seen, list.files(prep$target, pattern = "\\.asc$"))
  expect_true(all(out$diagnostics$output_validation$found == 4L))
  expect_s4_class(out$fcon, "SpatRaster")
})

test_that("invalid input manifests and geometry preserve outputs before Julia starts", {
  r <- make_test_raster(n = 3, crs = "EPSG:3857")
  prep <- conscape_prep(centersize = 2, buffer = 1, landmark = 1L,
                        r_target = r, r_mov = r, r_src = r, progress = FALSE)
  out <- tempfile("preserve-run-")
  dir.create(out)
  marker <- file.path(out, ".keep")
  writeLines("keep", marker)
  testthat::local_mocked_bindings(
    conscape_julia_start = function(...) stop("should not start Julia"),
    .package = "ConScapeRtools"
  )
  unlink(file.path(prep$src, "r_1.asc"))
  expect_error(run_conscape(conscape_prep = prep, out_dir = out, jl_home = "mock"),
               "identical sets")
  expect_true(file.exists(marker))
  expect_error(run_conscape(target_qualities = r, source_qualities = r,
                            affinities = "directory", out_dir = out, jl_home = "mock"),
               "single-layer")
  other <- r
  terra::crs(other) <- "EPSG:4326"
  expect_error(run_conscape(target_qualities = r, source_qualities = r,
                            affinities = other, out_dir = out, jl_home = "mock"), "share CRS")
  expect_true(file.exists(marker))
  expect_error(run_conscape(target_qualities = r, source_qualities = r,
                            affinities = r, out_dir = out, clear_dir = FALSE, jl_home = "mock"),
               "Files currently exist")
})

test_that("output paths cannot erase input data or broad directories", {
  r <- make_test_raster(n = 3, crs = "EPSG:3857")
  root <- tempfile("protected-")
  dir.create(root)
  source <- file.path(root, "source.tif")
  r <- terra::writeRaster(r, source)
  expect_error(conscape_prep(centersize = 2, buffer = 0, landmark = 1L,
                             asc_dir = root, clear_dir = TRUE,
                             r_target = r, r_mov = r, r_src = r), "overwrite input")
  expect_error(run_conscape(target_qualities = r, source_qualities = r, affinities = r,
                            out_dir = root, jl_home = "mock"), "overwrite input")
  expect_true(file.exists(source))
  for (path in c(getwd(), dirname(getwd()), tempdir(), path.expand("~"))) {
    expect_error(ConScapeRtools:::validate_conscape_output_dir(path, TRUE),
                 "dedicated working directory")
  }
})

test_that("Julia command arguments preserve spaces and code quotes", {
  received <- NULL
  testthat::local_mocked_bindings(system2 = function(command, args, ...) {
    received <<- args
    0L
  }, .package = "base")
  args <- c("--project=C:/analysis directory/Julia", "-e", 'println("two words")')
  ConScapeRtools:::run_julia_command("julia", args)
  expect_identical(received, shQuote(args))
})

test_that("Julia commands isolate Linux library paths without changing caller state", {
  old_path <- Sys.getenv("LD_LIBRARY_PATH", unset = NA_character_)
  on.exit(ConScapeRtools:::restore_conscape_envvar("LD_LIBRARY_PATH", old_path))
  Sys.setenv(LD_LIBRARY_PATH = "/R/library path:/other/library")
  received <- NULL
  testthat::local_mocked_bindings(
    Sys.info = function() c(sysname = platform),
    system2 = function(command, args, stdout, stderr, env) {
      received <<- list(command = command, args = args, env = env,
                        stdout = stdout, stderr = stderr)
      c("Julia output", "second line")
    }, .package = "base"
  )
  for (platform in c("Linux", "Windows", "Darwin")) {
    out <- ConScapeRtools:::run_julia_command("Julia bin/julia", "script path.jl",
                                            stdout = "output file", stderr = FALSE)
    expect_identical(out, c("Julia output", "second line"))
    expect_identical(received$command, "Julia bin/julia")
    expect_identical(received$args, shQuote("script path.jl"))
    expect_identical(received$env,
                     if (platform == "Linux") "LD_LIBRARY_PATH=''" else character())
    expect_identical(received$stdout, "output file")
    expect_false(received$stderr)
    expect_identical(Sys.getenv("LD_LIBRARY_PATH"), "/R/library path:/other/library")
  }
})

test_that("Julia command failures preserve status and the caller environment", {
  old_path <- Sys.getenv("LD_LIBRARY_PATH", unset = NA_character_)
  on.exit(ConScapeRtools:::restore_conscape_envvar("LD_LIBRARY_PATH", old_path))
  Sys.unsetenv("LD_LIBRARY_PATH")
  testthat::local_mocked_bindings(
    Sys.info = function() c(sysname = "Linux"),
    system2 = function(...) structure("Julia error", status = 1L), .package = "base"
  )
  out <- ConScapeRtools:::run_julia_command("julia")
  expect_identical(attr(out, "status"), 1L)
  expect_true(is.na(Sys.getenv("LD_LIBRARY_PATH", unset = NA_character_)))
  testthat::local_mocked_bindings(
    system2 = function(...) stop("launch failed"), .package = "base"
  )
  expect_error(ConScapeRtools:::run_julia_command("julia"), "launch failed")
  expect_true(is.na(Sys.getenv("LD_LIBRARY_PATH", unset = NA_character_)))
})

test_that("mosaics reject mixed geometry and multilayer tiles", {
  root <- tempfile("mixed-mosaic-")
  dir.create(root)
  x <- make_test_raster(n = 3, vals = 1, crs = "EPSG:3857")
  terra::writeRaster(x, file.path(root, "1.tif"))
  other <- x
  terra::crs(other) <- "EPSG:4326"
  terra::writeRaster(other, file.path(root, "2.tif"))
  expect_error(mosaic_conscape(root, tile_trim = 0), "matching CRS")
  terra::writeRaster(c(x, x), file.path(root, "2.tif"), overwrite = TRUE)
  expect_error(mosaic_conscape(root, tile_trim = 0), "single-layer")
})

test_that("prep reports effective alignment and partitions irregular target grids exactly", {
  r <- terra::rast(nrows = 7, ncols = 11, xmin = 0, xmax = 11,
                   ymin = 0, ymax = 7, vals = seq_len(77), crs = "EPSG:3857")
  terra::values(r)[c(1, 20, 77)] <- NA
  prep <- conscape_prep(centersize = 3.2, buffer = 0.2, landmark = 2L,
                        r_target = r, r_mov = r, r_src = r, progress = FALSE)
  expect_equal(prep$centersize, 4L)
  expect_equal(prep$buffer, 2L)
  expect_true(terra::same.crs(prep$cs_tiles, r))
  expect_true(terra::same.crs(terra::vect(file.path(prep$asc_dir, "tiles.shp")), r))
  merged <- mosaic_conscape(prep$target, tile_trim = prep$tile_trim, mask = !is.na(r), crs = terra::crs(r))
  expected <- terra::ifel(is.na(r), 0, r)
  expect_equal(terra::values(merged, mat = FALSE), terra::values(expected, mat = FALSE))
  prep2 <- conscape_prep(centersize = 4, buffer = 2, landmark = 2L,
                         r_target = r, r_mov = r, r_src = r, progress = FALSE)
  expect_false(identical(prep$asc_dir, prep2$asc_dir))
})

test_that("streamed mosaics retain precision, weighted means, coverage and backing files", {
  root <- tempfile("streamed-")
  dir.create(root)
  rasters <- lapply(1:7, function(i) {
    x <- make_test_raster(n = 3, vals = 1e8 + i / 10, crs = "EPSG:3857")
    v <- terra::values(x)
    v[1] <- NA
    if (i %% 2) v[2] <- NA
    terra::values(x) <- v
    terra::writeRaster(x, file.path(root, paste0(i, ".tif")), datatype = "FLT8S")
  })
  mat <- do.call(cbind, lapply(rasters, terra::values))
  for (size in c(1L, 2L, 3L, 64L)) {
    sum_out <- mosaic_conscape(root, tile_trim = 0, chunk_size = size)
    mean_out <- mosaic_conscape(root, tile_trim = 0, method = "mosaic", chunk_size = size)
    merge_out <- mosaic_conscape(root, tile_trim = 0, method = "merge", chunk_size = size)
    expect_equal(terra::values(sum_out, mat = FALSE), rowSums(mat, na.rm = TRUE), tolerance = 1e-14)
    expected <- rowMeans(mat, na.rm = TRUE)
    expected[is.nan(expected)] <- NA
    expect_equal(terra::values(mean_out, mat = FALSE), expected, tolerance = 1e-14)
    first <- apply(mat, 1, function(v) if (all(is.na(v))) NA_real_ else v[which(!is.na(v))[1]])
    expect_equal(terra::values(merge_out, mat = FALSE), first, tolerance = 1e-14)
    expect_true(all(file.exists(terra::sources(sum_out))))
    expect_equal(terra::values(terra::rast(terra::sources(sum_out))), terra::values(sum_out))
  }
})

test_that("criticality and exact calibration enforce their mathematical contracts", {
  r <- make_test_raster(n = 3, vals = 0.5)
  prep <- conscape_prep(centersize = 2, buffer = 1, landmark = 1L,
                        r_target = r, r_mov = r, r_src = r, progress = FALSE)
  expect_error(run_conscape(conscape_prep = prep, out_dir = tempfile(),
                            metrics = "criticality", jl_home = "mock"), "not additive")
  expect_error(run_conscape(target_qualities = r, source_qualities = r, affinities = r,
                            out_dir = tempfile(), metrics = "criticality", jl_home = "mock",
                            connectivity_function = "survival_probability"), "only supports")
  spec <- ConScapeRtools:::conscape_metric_specs()$criticality
  expect_identical(ConScapeRtools:::pick_mosaic_method(spec, "center"), "mosaic")
  for (cost in c(1e-7, 10, 1e6)) {
    scale <- ConScapeRtools:::calibrate_distance_scale(cost, 0.025)
    expect_equal(exp(-cost / scale), 0.025, tolerance = 1e-14)
  }
  expect_error(ConScapeRtools:::calibrate_distance_scale(0, 0.025), "finite positive")
})

test_that("scenario counts use classic padded geometry and effective alignment", {
  r <- make_test_raster(n = 11)
  scenarios <- window_design_scenarios(r, centersize = c(3, 6), buffer = 1, landmark = 2L)
  expect_equal(scenarios$centersize, c(6L, 4L))
  expect_equal(scenarios$buffer, c(2L, 2L))
  expect_equal(scenarios$expected_windows, c(4, 9))
  expect_error(window_design_scenarios(r, centersize = 3, buffer = 1, workers = 1.5), "positive integer")
  expect_error(conscape_sensitivity(one_out_of = 1.5), "positive integer")
  expect_error(conscape_sensitivity(target_equal_source = NA), "TRUE or FALSE")
})

test_that("automatic mosaicking refuses missing tile contributions", {
  mock_julia()
  r <- make_test_raster(n = 4)
  prep <- conscape_prep(centersize = 2, buffer = 1, landmark = 1L,
                        r_target = r, r_mov = r, r_src = r, progress = FALSE)
  testthat::local_mocked_bindings(juliaCall = function(...) "ok", .package = "ConScapeRtools")
  expect_error(run_conscape(conscape_prep = prep, out_dir = tempfile("missing-tiles-"),
                            parallel = TRUE, workers = 2L, jl_home = "mock"),
               "refusing to mosaic")
})

test_that("Julia startup does not create unnecessary module proxy finalizers", {
  testthat::local_mocked_bindings(
    check_julia_pkgs = function(...) TRUE,
    juliaImport = function(...) stop("unnecessary proxy import"),
    juliaEval = function(expr) {
      expect_match(expr, "using ConScape, SparseArrays, Statistics")
      invisible(NULL)
    }, .package = "ConScapeRtools"
  )
  expect_no_error(ConScapeRtools:::ConScapeR_setup("mock"))
})

test_that("efficiency proxies exclude unwritten tiles and reject fractional references", {
  r <- make_test_raster(n = 4)
  prep <- conscape_prep(centersize = 2, buffer = 1, landmark = 1L,
                        r_target = r, r_mov = r, r_src = r, progress = FALSE)
  diag <- prep$tile_diagnostics
  extra <- diag[1, , drop = FALSE]
  extra$written <- FALSE
  extra$n_positive_move <- extra$n_positive_target <- 1e6
  prep2 <- prep
  prep2$tile_diagnostics <- rbind(diag, extra)
  x <- conscape_efficiency_assessment(prep, prep2)
  expect_equal(x$estimated_source_target_cells[1], x$estimated_source_target_cells[2])
  for (reference in c(1.5, Inf, NaN, 1e20)) {
    expect_error(conscape_efficiency_assessment(prep, reference = reference),
                 "reference must identify")
  }
})
test_that("integer overflow is rejected before Julia starts", {
  r <- make_test_raster(n = 3, vals = 0.5)
  for (arg in c("blas_threads", "mosaic_chunk_size")) {
    args <- list(target_qualities = r, source_qualities = r, affinities = r,
                 out_dir = tempfile("overflow-"), jl_home = "unused")
    args[[arg]] <- .Machine$integer.max + 1
    expect_error(do.call(run_conscape, args), arg)
  }
})

test_that("workflow comparisons cannot hide missing coverage or failed outputs", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("examples", "validate_workflows.R", package = "ConScapeRtools"), env)
  r <- make_test_raster(n = 3, vals = 1)
  names(r) <- "fcon"
  expect_error(env$as_stack(NULL), "required workflow failed")
  missing <- r
  missing[1] <- NA
  expect_error(env$agreement(missing, r, "incomplete"), "coverage differs")
  blank <- r * NA_real_
  expect_error(env$agreement(blank, blank, "empty"), "No finite values")
  expect_true(env$agreement(r, r, "self")$identical_fp)
})
