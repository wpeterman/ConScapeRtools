## Parallel function -- JuliaConnectoR
# cs_par.func <- function(i,
#                         .jl_home,
#                         .hab_target,
#                         .hab_src,
#                         .mov_prob,
#                         .data_dir,
#                         .out_dir,
#                         .landmark,
#                         .theta,
#                         .exp_d){
#
#   iter <- paste0('-iter_', i)
#   r_target <- .hab_target[i]
#   r_source <- .hab_src[i]
#   r_res <- .mov_prob[i]
#
#   cs <- conscape(data_dir = .data_dir,
#                  out_dir = .out_dir,
#                  r_target = r_target,
#                  r_source = r_source,
#                  r_res = r_res,
#                  landmark = .landmark,
#                  theta = .theta,
#                  exp_d = .exp_d,
#                  iter = iter)
# }


# ConScape jl Function ----------------------------------------------------
#' ConScape jl function
#'
#' @description This is an internal function for running ConScape

#' @return A Julia function for running ConScape

#' @author Bill Peterman
#' @keywords internal

conscape <- function(data_dir,
                     out_dir,
                     r_target,
                     r_source,
                     r_res,
                     landmark,
                     theta,
                     exp_d,
                     iter){
  conscape <- JuliaConnectoR::juliaEval('function conscape(data_dir, out_dir, r_target, r_source, r_res, landmark, theta, exp_d, iter)

# set folders
datadir = joinpath(data_dir);
outdir_btwn = joinpath(out_dir, "btwn");
outdir_fcon = joinpath(out_dir, "fcon");

# created the output folders, if needed
if !isdir(out_dir)
    mkdir(out_dir)
end

if !isdir(outdir_btwn)
    mkdir(outdir_btwn)
end

if !isdir(outdir_fcon)
    mkdir(outdir_fcon)
end

# read habitat quality raster
hab_qual_target, meta_t = ConScape.readasc(joinpath(datadir, r_target));
hab_qual_source, meta_s = ConScape.readasc(joinpath(datadir, r_source));

## read movemement probability raster
mov_prob, meta_p = ConScape.readasc(joinpath(datadir, r_res));

keys(meta_p)


collect(values(meta_p))[1:end .!= 3]
collect(values(meta_p))[1:end .!= 3] == collect(values(meta_s))[1:end .!= 3]

non_matches = findall(xor.(isnan.(mov_prob), isnan.(hab_qual_target)))
mov_prob[non_matches] .= -9999
hab_qual_target[non_matches] .= -9999;

non_matches = findall(xor.(isnan.(mov_prob), isnan.(hab_qual_source)))
mov_prob[non_matches] .= -9999
hab_qual_source[non_matches] .= -9999;

adjacency_matrix = ConScape.graph_matrix_from_raster(mov_prob)
g = ConScape.Grid(size(mov_prob)...,
                    affinities = adjacency_matrix,
                    source_qualities = hab_qual_source,
                    target_qualities = ConScape.sparse(hab_qual_target),
                    costs = ConScape.mapnz(x -> -log(x), adjacency_matrix))

coarse_target_qualities = ConScape.coarse_graining(g, landmark)

g_coarse = ConScape.Grid(size(mov_prob)...,
    affinities=adjacency_matrix,
    source_qualities=hab_qual_source,
    target_qualities=coarse_target_qualities,
    costs=ConScape.mapnz(x -> -log(x), adjacency_matrix));
θ = theta

@time h_coarse = ConScape.GridRSP(g_coarse, θ=θ);

func_con = ConScape.connected_habitat(h_coarse,
        distance_transformation=x -> exp(-x/exp_d));


kbetw = ConScape.betweenness_kweighted(h_coarse,
                distance_transformation=x -> exp(-x/exp_d));

## Save results
ConScape.writeasc(joinpath(outdir_btwn, "betweenness" * iter * ".asc"), kbetw, meta_p);
ConScape.writeasc(joinpath(outdir_fcon, "fcon" * iter * ".asc"), func_con, meta_p);

end')

  conscape <- JuliaConnectoR::juliaFun("conscape")

  try(cs <- conscape(data_dir,
                     out_dir,
                     r_target,
                     r_source,
                     r_res,
                     landmark,
                     theta,
                     exp_d,
                     iter), silent = T)
}

## future availableCores plan
## future.apply future_lapply
## JuliaConnectoR juliaEval juliaFun juliaCall
