## R CMD check results

Checked locally on Windows 11 x64 with R 4.6.0 using:

```r
devtools::check(cran = TRUE, document = FALSE, vignettes = TRUE, manual = FALSE)
```

Status: OK (0 errors, 0 warnings, 0 notes).

On a first submission, CRAN's incoming checks may add a "New submission"
NOTE; that is expected.

## Release summary

This release continues the 0.4.x line that aligns `ConScapeRtools` with
the current ConScape API and stabilizes its tiled randomized
shortest-path workflow. Highlights since the last public release:

* Corrected `target_mode = "center"` so each tile's RSP target set
  really is the tile's interior cells (previously a no-op outside
  diagnostics). Per-tile contributions recombine via a new
  `mosaic_conscape(method = "sum")` reduction, reproducing the untiled
  ConScape solution to floating-point tolerance when the buffer covers
  the relevant landscape context, and converging exponentially in
  buffer otherwise.
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
* New integration tests
  (`test-tiled-untiled-agreement.R`,
  `test-tiled-untiled-convergence.R`,
  `test-sensitivity-convergence.R`) document buffer-truncation error
  and convergence behavior. These tests are gated behind
  `RUN_CONSCAPERTOOLS_INTEGRATION=true` and (for sensitivity)
  `RUN_CONSCAPERTOOLS_SENSITIVITY=true`, so they are skipped on CRAN.
* Renamed and rewrote the windowed/batch vignette to focus on
  validated workflows; the experimental ConScape dev backend is
  currently broken upstream and is now confined to a single
  "Experimental: Dev Backend (Currently Broken Upstream)" section
  with the captured error message and workaround.

## External software

ConScapeRtools interfaces with Julia and the Julia package ConScape for
its computational workflows. Unit tests, examples, and vignettes do not
launch Julia, so CRAN checks do not require an external Julia
installation or network access. Integration tests that do require Julia
are gated behind explicit environment variables and are skipped on
CRAN.
