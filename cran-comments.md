## R CMD check results

Checked locally on Windows 11 x64 using R 4.5.3 with:

```r
R CMD check --as-cran --no-manual ConScapeRtools_0.4.15.tar.gz
```

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* New submission.

This is expected because ConScapeRtools is not currently on CRAN.

## Release summary

This release refreshes ConScapeRtools for the current ConScape API and tiling
workflow. User-facing changes include:

* ConScape-aligned argument names for target qualities, source qualities,
  affinities, cost functions, connectivity functions, metrics, and sensitivity
  outputs, while preserving legacy aliases where possible.
* Sensitivity specification helpers and Julia runner support for ConScape
  sensitivity outputs when they are available in the installed Julia package.
* Promotion of `distance_scale` as the user-facing exponential distance decay
  parameter, with `exp_d` retained for backwards compatibility.
* More explicit `tile_design()` diagnostics for calibrated decay, tile width,
  tile trim, expected tile counts, overlap area factors, and trim proximity.

## External software

ConScapeRtools interfaces with Julia and the Julia package ConScape for
computational workflows. Unit tests and examples avoid launching Julia so CRAN
checks do not require external Julia installation or network access.
