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
