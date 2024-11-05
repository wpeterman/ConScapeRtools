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

mov_tile <- tile_rast(r = resist,
                      make_tiles = target_tile ,
                      out_dir = file.path(target_tile$asc_dir, 'mov'),
                      clear_dir = T)
