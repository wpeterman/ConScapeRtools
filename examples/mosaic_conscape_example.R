library(ConScapeRtools)

## Import data
s <- system.file("data/suitability.asc", package = "ConScapeRtools")
source <- terra::rast(s)

a <- system.file("data/affinity.asc", package = "ConScapeRtools")
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
landmark <- 2L # Must be an integer, not numeric

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
cs_run <- run_conscape(out_dir = file.path(prep$asc_dir,"results"),
                       conscape_prep = prep,
                       theta = theta,
                       exp_d = exp_d,
                       jl_home = jl_home,
                       parallel = T,
                       workers = 5,
                       end_multisession = T)

## Put output tiles together
cs_btwn <- mosaic_conscape(asc_dir = cs_run$outdir_btwn,
                           tile_trim = tile_trim)
cs_fcon <- mosaic_conscape(asc_dir = cs_run$outdir_fcon,
                           tile_trim = tile_trim)
plot(c(cs_btwn, cs_fcon))
