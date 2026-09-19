## R CMD check results

Checked from the built source archive on Windows 11 x64 with:

* R 4.6.1, `R CMD check --as-cran`: 0 errors, 0 warnings, 1 note.
* R 4.5.3, `R CMD check --as-cran --no-manual`: status OK.

The R 4.6.1 note is the expected `New submission` note. The PDF and HTML
manuals, examples, tests, vignettes, and vignette rebuilding all passed in the
full R 4.6.1 check.

## Release summary

This release prepares the package for its first CRAN submission. Highlights
include:

* Corrected `target_mode = "center"` so each tile's RSP target set
  really is the tile's interior cells (previously a no-op outside
  diagnostics). Per-tile contributions recombine via a new
  `mosaic_conscape(method = "sum")` reduction, reproducing the untiled
  ConScape solution to floating-point tolerance when the buffer covers
  the relevant landscape context. Empirical integration tests verify that
  error decreases as buffer size increases.
* Per-output-layer mosaic dispatch: additive metrics (`fcon`, `btwn`,
  `btwn_qweighted`, `criticality`) use sum mosaic under center mode and
  mean mosaic under legacy full mode; sensitivity surfaces are
  tile-local landscape-summary derivatives and always use mean mosaic.
* `run_conscape()` refuses sensitivity requests against
  center-mode preps, because center mode violates ConScape's
  `target_equal_source = TRUE` precondition inside every tile. The
  refusal happens before Julia starts and points the user to the
  documented workaround (`target_mode = "full"` with `landmark = 1L`).
* New `inst/examples/validate_workflows.R` harness reports identity
  and timing across untiled, classic tiled, center-buffer, and legacy
  full-target runs. Its output is embedded in the new vignette
  `Tiled ConScape: Validation and Performance`.
* Integration tests
  (`test-tiled-untiled-agreement.R`,
  `test-tiled-untiled-convergence.R`,
  `test-sensitivity-convergence.R`) document buffer-truncation error
  and convergence behavior. The installed Julia 1.12.6 and ConScape 0.3.0
  configuration passed 40 integration assertions for the threaded,
  distributed, agreement, and convergence workflows. These tests are gated behind
  `RUN_CONSCAPERTOOLS_INTEGRATION=true` and (for sensitivity)
  `RUN_CONSCAPERTOOLS_SENSITIVITY=true`, so they are skipped on CRAN.
* Renamed and rewrote the windowed/batch vignette to focus on
  validated workflows; the experimental ConScape dev backend is
  currently broken upstream and is now confined to a single
  "Experimental: Dev Backend (Currently Broken Upstream)" section
  with the captured error message and workaround.
* Removed private `JuliaConnectoR` state access, made Julia package installation
  explicitly opt-in, restored caller environment and parallel state, and added
  strict whole-number validation for discrete arguments.
* Added runnable examples, package-level documentation, spelling and URL
  validation, direct development-backend coverage, and a cross-platform check
  workflow for Windows, macOS, R release, R devel, and R oldrel-1.

## External software

ConScapeRtools interfaces with Julia and the Julia package ConScape for
its computational workflows. Unit tests, examples, and vignettes do not
launch Julia, so CRAN checks do not require an external Julia
installation or network access. Integration tests that require Julia
are gated behind explicit environment variables and are skipped on
CRAN.
