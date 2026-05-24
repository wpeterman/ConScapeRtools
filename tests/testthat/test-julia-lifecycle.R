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
