library(ConScapeRtools)

## Import data
s <- system.file("data/suitability.asc", package = "ConScapeRtools")
source <- terra::rast(s)

a <- system.file("data/affinity.asc", package = "ConScapeRtools")
resist <- terra::rast(a)

jl_home <- "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.10.5/bin/"

## Run tile design to get parameters
td <- tile_design(r_mov = resist,
                  r_target = source,
                  max_d = 7000,
                  theta = 0.1,
                  jl_home = jl_home)

target_tile <- make_tiles(tile_d = td$tile_d,
                          tile_trim = td$tile_trim,
                          r = source,
                          clear_dir = T)
