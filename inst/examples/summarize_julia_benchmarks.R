# Rscript summarize_julia_benchmarks.R BENCHMARK_LOG_DIRECTORY OUTPUT_CSV
# Summarize raw logs from benchmark_julia_windows.ps1 without editing numbers.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, dir.exists(args[1]), !file.exists(args[2]))
peaks <- read.csv(file.path(args[1], "peaks.csv"))
rows <- lapply(seq_len(nrow(peaks)), function(i) {
  x <- read.csv(file.path(args[1], paste0(peaks$Case[i], ".csv")), comment.char = "#")
  stopifnot(nrow(x) == 3L, all(is.finite(x$seconds)), all(x$seconds > 0))
  data.frame(method = x$method[1], side = x$side[1], targets = x$targets[1],
              blas = x$blas[1], block = x$block[1],
              seconds_median = median(x$seconds),
              allocated_mb_median = median(x$allocated_mb),
              peak_working_set_mb = peaks$PeakWorkingSetMB[i])
})
write.csv(do.call(rbind, rows), args[2], row.names = FALSE)
