# Run from a source checkout after devtools::load_all(). No Julia runs at build time.
# Rscript inst/examples/refresh_validation.R OUTPUT_DIRECTORY
# Set CONSCAPERTOOLS_JL_HOME first. OUTPUT_DIRECTORY must not already exist and
# must be separate from all disposable prep/run directories.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, !dir.exists(args[1]))
dir.create(args[1], recursive = TRUE)
output <- normalizePath(args[1], winslash = "/")
pkgload::load_all()
source(system.file("examples", "validate_workflows.R", package = "ConScapeRtools"))
habitat <- terra::rast(system.file("extdata", "suitability.asc", package = "ConScapeRtools"))
affinity <- terra::rast(system.file("extdata", "affinity.asc", package = "ConScapeRtools"))
td <- tile_design(r_mov = affinity, r_target = habitat, max_d = 1200,
                  theta = 0.1, jl_home = Sys.getenv("CONSCAPERTOOLS_JL_HOME"),
                  trim_threshold = 0.05, overlap_area_factor = 2, landmark = 5L)
saveRDS(td, file.path(output, "calibration.rds"))
capture.output(summary(td), file = file.path(output, "calibration.txt"))
for (buffer in unique(c(td$buffer, 20L))) {
  result <- validate_workflows(centersize = td$centersize, buffer = buffer,
                               theta = td$theta, distance_scale = td$distance_scale,
                               landmark = td$landmark)
  write.csv(result$summary, file.path(output, paste0("agreement-buffer-", buffer, ".csv")),
            row.names = FALSE)
  write.csv(result$timings, file.path(output, paste0("timing-buffer-", buffer, ".csv")),
            row.names = FALSE)
}
conscape_julia_start(Sys.getenv("CONSCAPERTOOLS_JL_HOME"), quiet = TRUE)
backend <- JuliaConnectoR::juliaEval(paste0(
  "using Pkg; p = only(p for p in values(Pkg.dependencies()) if p.name == \"ConScape\"); ",
  "string(VERSION, \"; ConScape \", p.version, \"; tree \", p.tree_hash)"))
conscape_julia_stop()
writeLines(c(paste("Recorded:", Sys.time()), paste("Julia:", backend),
             capture.output(sessionInfo()), "Input and source MD5:",
             capture.output(tools::md5sum(c("inst/extdata/conscape.jl",
               "inst/extdata/suitability.asc", "inst/extdata/affinity.asc")))),
           file.path(output, "provenance.txt"))
