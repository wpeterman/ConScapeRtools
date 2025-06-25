library(ConScapeRtools)

## Import data
s <- system.file("extdata", "suitability.asc", package = "ConScapeRtools")
source <- terra::rast(s)

a <- system.file("extdata", "affinity.asc", package = "ConScapeRtools")
resist <- terra::rast(a)

jl_home <- "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.10.5/bin/"

td <- tile_design(r_mov = resist,
                  r_target = source,
                  max_d = 7000,
                  theta = 0.1,
                  jl_home = jl_home)

## Tile dimension
tile_d <- td$tile_d

# How much to trim tiles
tile_trim <- td$tile_trim

# Makes computation more efficient
landmark <- 5L # Must be an integer, not numeric

# Controls level of randomness of paths
theta <- td$theta

# Controls rate of decay with distance
exp_d <- td$exp_d

## Prepare data for analysis
prep <- conscape_prep(tile_d = tile_d,
                      tile_trim = tile_trim,
                      r_target = source,
                      r_mov = resist,
                      r_src = source,
                      clear_dir = T,
                      landmark = landmark)

## Run ConScape
## No parallelization
cs_run.serial <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
                              conscape_prep = prep,
                              theta = theta,
                              exp_d = exp_d,
                              jl_home = jl_home,
                              parallel = F)

## Threaded parallel
cs_run.thread <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
                            conscape_prep = prep,
                            theta = theta,
                            exp_d = exp_d,
                            jl_home = jl_home,
                            parallel = T,
                            workers = 4)

## Threaded parallel
cs_run.dist <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
                           conscape_prep = prep,
                           theta = theta,
                           exp_d = exp_d,
                           jl_home = jl_home,
                           parallel = T,
                           workers = 4,
                           distributed = TRUE)
plot(cs_run.dist)

## No tiling --> Only attempt with small to moderate sized rasters!
cs_run <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
                       hab_target = source,
                       hab_src = source,
                       mov_prob = resist,
                       theta = theta,
                       exp_d = exp_d,
                       jl_home = jl_home)

