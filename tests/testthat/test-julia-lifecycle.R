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
