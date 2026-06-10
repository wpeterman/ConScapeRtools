test_that("JuliaConnectoR finalized refs are cleared before startup", {
  ns <- asNamespace("JuliaConnectoR")
  pkg_local <- get("pkgLocal", envir = ns, inherits = FALSE)
  old_refs <- pkg_local$finalizedRefs
  on.exit(pkg_local$finalizedRefs <- old_refs, add = TRUE)

  pkg_local$finalizedRefs <- as.raw(seq_len(8))
  cleared <- getFromNamespace(
    "clear_juliaconnector_finalized_refs",
    "ConScapeRtools"
  )()

  expect_true(cleared)
  expect_null(pkg_local$finalizedRefs)
})

test_that("stop_conscape_julia clears stale finalized refs after stopping", {
  ns <- asNamespace("JuliaConnectoR")
  pkg_local <- get("pkgLocal", envir = ns, inherits = FALSE)
  old_refs <- pkg_local$finalizedRefs
  on.exit(pkg_local$finalizedRefs <- old_refs, add = TRUE)

  called <- FALSE
  testthat::local_mocked_bindings(
    stopJulia = function() {
      called <<- TRUE
      invisible(TRUE)
    },
    .env = asNamespace("ConScapeRtools")
  )

  pkg_local$finalizedRefs <- as.raw(seq_len(8))
  getFromNamespace("stop_conscape_julia", "ConScapeRtools")()

  expect_true(called)
  expect_null(pkg_local$finalizedRefs)
})

test_that("package Julia calls clear stale finalized refs", {
  ns <- asNamespace("JuliaConnectoR")
  pkg_local <- get("pkgLocal", envir = ns, inherits = FALSE)
  old_refs <- pkg_local$finalizedRefs
  on.exit(pkg_local$finalizedRefs <- old_refs, add = TRUE)

  saw_cleared_refs <- FALSE
  testthat::local_mocked_bindings(
    juliaCall = function(...) {
      saw_cleared_refs <<- is.null(pkg_local$finalizedRefs)
      "ok"
    },
    .env = asNamespace("ConScapeRtools")
  )

  pkg_local$finalizedRefs <- as.raw(seq_len(8))
  out <- getFromNamespace("juliaCall_conscape", "ConScapeRtools")("identity", 1)

  expect_equal(out, "ok")
  expect_true(saw_cleared_refs)
})

test_that("package Julia calls clear stale finalized refs after errors", {
  ns <- asNamespace("JuliaConnectoR")
  pkg_local <- get("pkgLocal", envir = ns, inherits = FALSE)
  old_refs <- pkg_local$finalizedRefs
  on.exit(pkg_local$finalizedRefs <- old_refs, add = TRUE)

  testthat::local_mocked_bindings(
    juliaCall = function(...) {
      pkg_local$finalizedRefs <<- as.raw(seq_len(8))
      stop("Julia call failed", call. = FALSE)
    },
    .env = asNamespace("ConScapeRtools")
  )

  expect_error(
    getFromNamespace("juliaCall_conscape", "ConScapeRtools")("identity", 1),
    "Julia call failed"
  )
  expect_null(pkg_local$finalizedRefs)
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
    .env = asNamespace("ConScapeRtools")
  )

  expect_invisible(conscape_julia_start("C:/Julia/bin", quiet = TRUE))
  expect_true(started)

  expect_invisible(conscape_julia_stop())
  expect_true(stopped)
})

test_that("conscape_julia_status does not start Julia", {
  setup_checked <- FALSE
  testthat::local_mocked_bindings(
    juliaSetupOk = function() {
      setup_checked <<- TRUE
      TRUE
    },
    .env = asNamespace("ConScapeRtools")
  )

  expect_false(conscape_julia_status())
  expect_false(setup_checked)
})

test_that("conscape_julia_start reports invalid setup", {
  testthat::local_mocked_bindings(
    juliaSetupOk = function() FALSE,
    .env = asNamespace("ConScapeRtools")
  )

  expect_error(
    conscape_julia_start("C:/Julia/bin", quiet = TRUE),
    "path to the Julia"
  )
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
