test_that("stop_conscape_julia uses the public JuliaConnectoR API", {
  called <- FALSE
  testthat::local_mocked_bindings(
    stopJulia = function() {
      called <<- TRUE
      invisible(TRUE)
    },
    .package = "ConScapeRtools"
  )

  getFromNamespace("stop_conscape_julia", "ConScapeRtools")()

  expect_true(called)
  expect_false(conscape_julia_status())
})

test_that("package Julia calls use the public juliaCall wrapper", {
  seen <- NULL
  testthat::local_mocked_bindings(
    juliaCall = function(...) {
      seen <<- list(...)
      "ok"
    },
    .package = "ConScapeRtools"
  )

  out <- getFromNamespace("juliaCall_conscape", "ConScapeRtools")("identity", 1)

  expect_equal(out, "ok")
  expect_equal(seen, list("identity", 1))
})

test_that("package Julia calls preserve JuliaConnectoR errors", {
  testthat::local_mocked_bindings(
    juliaCall = function(...) stop("Julia call failed", call. = FALSE),
    .package = "ConScapeRtools"
  )

  expect_error(
    getFromNamespace("juliaCall_conscape", "ConScapeRtools")("identity", 1),
    "Julia call failed"
  )
})

test_that("conscape_julia_status does not start Julia", {
  setup_checked <- FALSE
  testthat::local_mocked_bindings(
    stopJulia = function() invisible(TRUE),
    juliaSetupOk = function() {
      setup_checked <<- TRUE
      TRUE
    },
    .package = "ConScapeRtools"
  )

  conscape_julia_stop()
  expect_false(conscape_julia_status())
  expect_false(setup_checked)
})

test_that("conscape_julia_start reports invalid setup", {
  testthat::local_mocked_bindings(
    juliaSetupOk = function() FALSE,
    .package = "ConScapeRtools"
  )

  expect_error(
    conscape_julia_start("C:/Julia/bin", quiet = TRUE),
    "path to the Julia"
  )
})

test_that("Julia setup never installs packages unless explicitly requested", {
  expressions <- character()
  testthat::local_mocked_bindings(
    juliaEval = function(expr, ...) {
      expressions <<- c(expressions, expr)
      FALSE
    },
    .package = "ConScapeRtools"
  )

  expect_error(
    getFromNamespace("ConScapeR_setup", "ConScapeRtools")(
      "C:/Julia/bin",
      install_libraries = FALSE
    ),
    "not installed"
  )
  expect_false(any(grepl("Pkg.add", expressions, fixed = TRUE)))
})

test_that("explicit Julia installation adds ConScape and restores JULIA_BINDIR", {
  expressions <- character()
  old_bindir <- Sys.getenv("JULIA_BINDIR", unset = NA_character_)
  testthat::local_mocked_bindings(
    juliaEval = function(expr, ...) {
      expressions <<- c(expressions, expr)
      TRUE
    },
    juliaImport = function(pkg) pkg,
    .package = "ConScapeRtools"
  )

  expect_invisible(
    getFromNamespace("ConScapeR_setup", "ConScapeRtools")(
      "C:/Julia/bin",
      install_libraries = TRUE
    )
  )
  expect_equal(sum(grepl("Pkg.add", expressions, fixed = TRUE)), 1L)
  expect_match(expressions[grepl("Pkg.add", expressions, fixed = TRUE)], "ConScape")
  expect_identical(Sys.getenv("JULIA_BINDIR", unset = NA_character_), old_bindir)
})

test_that("exported Julia helpers start and stop", {
  started <- FALSE
  stopped <- FALSE

  testthat::local_mocked_bindings(
    juliaSetupOk = function() TRUE,
    ConScapeR_setup = function(julia_path, install_libraries = FALSE) {
      started <<- TRUE
      expect_equal(julia_path, "C:/Julia/bin")
      expect_false(install_libraries)
      invisible(TRUE)
    },
    stopJulia = function() {
      stopped <<- TRUE
      invisible(TRUE)
    },
    .package = "ConScapeRtools"
  )

  old_bindir <- Sys.getenv("JULIA_BINDIR", unset = NA_character_)
  expect_invisible(conscape_julia_start("C:/Julia/bin", quiet = TRUE))
  expect_true(started)
  expect_true(conscape_julia_status())
  expect_identical(Sys.getenv("JULIA_BINDIR", unset = NA_character_), old_bindir)

  expect_invisible(conscape_julia_stop())
  expect_true(stopped)
})

test_that("conscape_dev_backend_setup prepares an isolated Julia project", {
  julia_bin <- file.path(tempdir(), "fake-julia-bin")
  dir.create(julia_bin, recursive = TRUE, showWarnings = FALSE)
  julia_exe <- file.path(julia_bin, "julia.exe")
  writeLines("", julia_exe)

  project <- file.path(tempdir(), "conscape-dev-project-test")
  calls <- list()
  testthat::local_mocked_bindings(
    run_julia_command = function(command, args, stdout = TRUE, stderr = TRUE) {
      calls[[length(calls) + 1L]] <<- list(
        command = command,
        args = args,
        stdout = stdout,
        stderr = stderr
      )
      "ConScape dev backend project prepared"
    },
    .env = asNamespace("ConScapeRtools")
  )

  out <- conscape_dev_backend_setup(
    jl_home = julia_bin,
    project = project,
    rev = "alg_efficiency",
    url = "https://github.com/ConScape/ConScape.jl",
    quiet = TRUE
  )

  expect_equal(out, normalizePath(project, winslash = "/", mustWork = TRUE))
  expect_equal(calls[[1]]$command, julia_exe)
  expect_true("--startup-file=no" %in% calls[[1]]$args)
  expect_true("alg_efficiency" %in% calls[[1]]$args)
  expect_true("https://github.com/ConScape/ConScape.jl" %in% calls[[1]]$args)
})

test_that("conscape_dev_backend_setup reports Julia setup failures", {
  julia_bin <- file.path(tempdir(), "fake-julia-bin-fail")
  dir.create(julia_bin, recursive = TRUE, showWarnings = FALSE)
  writeLines("", file.path(julia_bin, "julia.exe"))

  testthat::local_mocked_bindings(
    run_julia_command = function(...) {
      out <- "Pkg failed"
      attr(out, "status") <- 1L
      out
    },
    .env = asNamespace("ConScapeRtools")
  )

  expect_error(
    conscape_dev_backend_setup(
      jl_home = julia_bin,
      project = file.path(tempdir(), "conscape-dev-project-fail"),
      quiet = TRUE
    ),
    "Failed to prepare"
  )
})

test_that("conscape_sensitivity_setup installs the sensitivity branch and parses the version", {
  julia_bin <- file.path(tempdir(), "fake-julia-bin-sens")
  dir.create(julia_bin, recursive = TRUE, showWarnings = FALSE)
  julia_exe <- file.path(julia_bin, "julia.exe")
  writeLines("", julia_exe)

  calls <- list()
  testthat::local_mocked_bindings(
    run_julia_command = function(command, args, stdout = TRUE, stderr = TRUE) {
      calls[[length(calls) + 1L]] <<- list(command = command, args = args)
      c("Resolving package versions...",
        "ConScape sensitivity version: 0.3.0")
    },
    .env = asNamespace("ConScapeRtools")
  )

  version <- conscape_sensitivity_setup(jl_home = julia_bin, quiet = TRUE)

  expect_identical(version, "0.3.0")
  expect_equal(calls[[1]]$command, julia_exe)
  expect_true("--startup-file=no" %in% calls[[1]]$args)
  expect_true("sensitivity" %in% calls[[1]]$args)
  expect_true("https://github.com/ConScape/ConScape.jl" %in% calls[[1]]$args)
  # force defaults to FALSE -> passes "false" as the final positional arg
  expect_true("false" %in% calls[[1]]$args)
  expect_false("true" %in% calls[[1]]$args)
})

test_that("conscape_sensitivity_setup passes force = TRUE through to Julia", {
  julia_bin <- file.path(tempdir(), "fake-julia-bin-sens-force")
  dir.create(julia_bin, recursive = TRUE, showWarnings = FALSE)
  writeLines("", file.path(julia_bin, "julia.exe"))

  calls <- list()
  testthat::local_mocked_bindings(
    run_julia_command = function(command, args, stdout = TRUE, stderr = TRUE) {
      calls[[length(calls) + 1L]] <<- list(args = args)
      "ConScape sensitivity version: 0.3.0"
    },
    .env = asNamespace("ConScapeRtools")
  )

  conscape_sensitivity_setup(jl_home = julia_bin, force = TRUE, quiet = TRUE)
  expect_true("true" %in% calls[[1]]$args)
})

test_that("conscape_sensitivity_setup reports Julia install failures", {
  julia_bin <- file.path(tempdir(), "fake-julia-bin-sens-fail")
  dir.create(julia_bin, recursive = TRUE, showWarnings = FALSE)
  writeLines("", file.path(julia_bin, "julia.exe"))

  testthat::local_mocked_bindings(
    run_julia_command = function(...) {
      out <- "Installed ConScape is missing sensitivity API symbols: sensitivity"
      attr(out, "status") <- 1L
      out
    },
    .env = asNamespace("ConScapeRtools")
  )

  expect_error(
    conscape_sensitivity_setup(jl_home = julia_bin, quiet = TRUE),
    "Failed to install the ConScape sensitivity build"
  )
})

test_that("conscape_sensitivity_setup validates its arguments", {
  expect_error(conscape_sensitivity_setup(jl_home = ""),
               "jl_home must be a single non-empty")
  expect_error(conscape_sensitivity_setup(jl_home = "x", rev = ""),
               "rev must be a single non-empty")
  expect_error(conscape_sensitivity_setup(jl_home = "x", url = NA_character_),
               "url must be a single non-empty")
  expect_error(conscape_sensitivity_setup(jl_home = "x", force = "yes"),
               "force must be TRUE or FALSE")
})
