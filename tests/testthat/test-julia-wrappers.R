test_that("mat2rast preserves matrix values on a raster template", {
  mat2rast <- getFromNamespace("mat2rast", "ConScapeRtools")
  template <- make_test_raster(n = 2, vals = 0, crs = "EPSG:4326")
  mat <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)

  out <- mat2rast(mat, template)

  expect_s4_class(out, "SpatRaster")
  expect_equal(as.vector(terra::ext(out)), as.vector(terra::ext(template)))
  expect_equal(terra::crs(out), terra::crs(template))
  expect_equal(as.vector(terra::as.matrix(out, wide = TRUE)), as.vector(mat))
})

test_that("Grid prepares raster inputs and dispatches to Julia", {
  Grid <- getFromNamespace("Grid", "ConScapeRtools")
  calls <- list()
  testthat::local_mocked_bindings(
    juliaLet = function(code, ...) {
      calls[[length(calls) + 1L]] <<- list(code = code, args = list(...))
      structure(list(code = code), class = "julia_grid")
    },
    .env = asNamespace("ConScapeRtools")
  )

  r <- make_test_raster(n = 2, vals = c(1, NaN, 3, 4))
  out <- Grid(affinities = r, sources = r, targets = r, costs = "x -> -log(x)")

  expect_s3_class(out, "julia_grid")
  expect_match(calls[[1]]$code, "MinusLog")
  expect_true(is.matrix(calls[[1]]$args$affinities))
  expect_true(is.nan(calls[[1]]$args$sources[1, 2]))
})

test_that("Grid handles matrix costs separately from transformed affinities", {
  Grid <- getFromNamespace("Grid", "ConScapeRtools")
  captured <- NULL
  testthat::local_mocked_bindings(
    juliaLet = function(code, ...) {
      captured <<- list(code = code, args = list(...))
      "grid"
    },
    .env = asNamespace("ConScapeRtools")
  )

  affinities <- matrix(c(1, NaN, 3, 4), 2)
  sources <- matrix(1, 2, 2)
  targets <- matrix(1, 2, 2)
  costs <- matrix(2, 2, 2)

  expect_equal(Grid(affinities, sources, targets, costs), "grid")
  expect_match(captured$code, "graph_matrix_from_raster\\(costs\\)")
  expect_true(is.nan(captured$args$costs[2, 1]))
  expect_true(is.nan(captured$args$sources[2, 1]))
})

test_that("Julia wrapper expressions suppress informational logging", {
  wrap_expr <- getFromNamespace("julia_warn_or_error_expr", "ConScapeRtools")
  out <- wrap_expr("ConScape.GridRSP(g, theta = theta)")

  expect_match(out, "Logging.ConsoleLogger\\(stderr, Logging.Warn\\)")
  expect_match(out, "ConScape.GridRSP")
})

test_that("Julia wrapper functions pass expected code to juliaLet", {
  wrappers <- list(
    vec2mat = list(fun = getFromNamespace("vec2mat", "ConScapeRtools"), args = list(vec = 1:3, g = "g"), pattern = "_vec_to_grid"),
    GridRSP = list(fun = getFromNamespace("GridRSP", "ConScapeRtools"), args = list(g = "g", theta = 0.1), pattern = "GridRSP"),
    betweenness_qweighted = list(fun = getFromNamespace("betweenness_qweighted", "ConScapeRtools"), args = list(h = "h"), pattern = "betweenness_qweighted"),
    betweenness_kweighted = list(fun = getFromNamespace("betweenness_kweighted", "ConScapeRtools"), args = list(h = "h", alpha = 2), pattern = "betweenness_kweighted"),
    coarse_grid = list(fun = getFromNamespace("coarse_grid", "ConScapeRtools"), args = list(g = "g", land_mark = 5L), pattern = "coarse_graining"),
    connected_habitat = list(fun = getFromNamespace("connected_habitat", "ConScapeRtools"), args = list(h = "h", alpha = 2), pattern = "connected_habitat"),
    expected_cost = list(fun = getFromNamespace("expected_cost", "ConScapeRtools"), args = list(h = "h"), pattern = "expected_cost"),
    sensitivity = list(fun = getFromNamespace("sensitivity", "ConScapeRtools"), args = list(h = "h", wrt = "Q"), pattern = "sensitivity")
  )

  code_seen <- character()
  testthat::local_mocked_bindings(
    juliaLet = function(code, ...) {
      code_seen <<- c(code_seen, code)
      paste0("ok-", length(code_seen))
    },
    .env = asNamespace("ConScapeRtools")
  )

  for (wrapper in wrappers) {
    result <- do.call(wrapper$fun, wrapper$args)
    expect_match(result, "^ok-")
    expect_match(tail(code_seen, 1), wrapper$pattern)
  }
})

