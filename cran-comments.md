## R CMD check results

ConScapeRtools 0.4.32 was checked from its built source archive on Windows 11
x64 with R 4.6.1 and `R CMD check --as-cran`: 0 errors, 0 warnings, 1 note.
The note is the expected `New submission` note. The full check includes
the PDF and HTML manuals, examples, tests, and vignette rebuilding.

The validation-branch workflow additionally checks Windows and macOS R-release,
plus Linux R-release, R-devel, and R-oldrel-1. These checks must pass before
merging the branch or preparing a submission. Real Julia integration tests run
separately on the three R-release platforms.

## Release summary

Version 0.4.32 declares GPL-2 for the combined package and credits the
ConScape.jl contributors for the adapted numerical routines. The installed
COPYRIGHTS file records the upstream files and versions, modifications, and
the retained MIT notice for preexisting code. COPYING supplies the GPL-2 text
in the source distribution; installed users can use the copy supplied with R.

This release prepares the supported classic backend for its first CRAN
submission:

* Default expected-cost connectivity and betweenness use bounded target blocks
  and one sparse factorization per tile. Numerical comparisons against upstream
  ConScape cover masked landscapes, dense and sparse targets, unequal source
  and target quality, cost transformations, theta values, and block sizes.
* Tile preparation avoids full-landscape extension. Disk-backed mosaics preserve
  sum and weighted-mean semantics across bounded groups of tiles. Inputs are
  validated before clearing output directories, and incomplete runs cannot
  silently produce mosaics.
* Additive center-target metrics use sum mosaics; finite buffers remain an
  approximation. Full-target metrics and sensitivity surfaces use trimmed mean
  mosaics. Criticality is not additive over target partitions, so center-target
  criticality is rejected. Sensitivity requires full targets and landmark one.
* The validation vignette includes regenerated calibration and tiled-versus-
  untiled comparisons, reproducible benchmark scripts, raw measurements,
  environment records, and explicit limits on performance and accuracy claims.
* The experimental development backend is disabled before Julia startup or
  output deletion. Its arguments remain available for a future repair, and the
  vignette records the failing upstream revision and requirements for reenabling
  it. The classic backend remains available.
* Julia dependency installation is opt-in. Runs restore caller state and remove
  only workers they created.

## External software and tests

ConScapeRtools interfaces with Julia and the Julia package ConScape. CRAN's unit
tests, runnable examples, and vignette builds do not launch Julia, install
external software, or require network access. Live integration tests require
explicit `RUN_CONSCAPERTOOLS_INTEGRATION=true`; sensitivity tests additionally
require `RUN_CONSCAPERTOOLS_SENSITIVITY=true` and the documented upstream
sensitivity build. Both groups are skipped on CRAN.

Local numerical validation used Julia 1.12.6 with registered ConScape 0.3.0
and the documented sensitivity build. The 304 direct Julia assertions passed
with each build. Performance measurements are workload-specific, not a claim
of optimal performance on every landscape or computer.

This is a first CRAN submission; there are no CRAN reverse dependencies.
