## R CMD check results

Checked locally on Windows 11 x64 using R 4.5.3 with:

```r
R CMD check --as-cran --no-manual ConScapeRtools_0.4.7.tar.gz
```

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* New submission.

This is expected because ConScapeRtools is not currently on CRAN.

## External software

ConScapeRtools interfaces with Julia and the Julia package ConScape for
computational workflows. Unit tests and examples avoid launching Julia so CRAN
checks do not require external Julia installation or network access.
