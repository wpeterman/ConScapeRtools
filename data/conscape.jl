function conscape(src_dir, mov_dir, target_dir, out_dir, r_target, r_source, r_res, land_mark, theta, exp_d, NA_val, iter)

# set folders
# datadir = joinpath(target_dir);
targetdir = joinpath(target_dir);
movdir = joinpath(mov_dir);
srcdir = joinpath(src_dir);
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
hab_qual_target, meta_t = ConScape.readasc(joinpath(targetdir, r_target));
hab_qual_source, meta_s = ConScape.readasc(joinpath(srcdir, r_source));

## read movemement probability raster
mov_prob, meta_p = ConScape.readasc(joinpath(movdir, r_res));

keys(meta_p)


collect(values(meta_p))[1:end .!= 3]
collect(values(meta_p))[1:end .!= 3] == collect(values(meta_s))[1:end .!= 3]

non_matches = findall(xor.(isnan.(mov_prob), isnan.(hab_qual_target)))
mov_prob[non_matches] .= NA_val
hab_qual_target[non_matches] .= NA_val;

non_matches = findall(xor.(isnan.(mov_prob), isnan.(hab_qual_source)))
mov_prob[non_matches] .= NA_val
hab_qual_source[non_matches] .= NA_val;

adjacency_matrix = ConScape.graph_matrix_from_raster(mov_prob)
g = ConScape.Grid(size(mov_prob)...,
                    affinities = adjacency_matrix,
                    source_qualities = hab_qual_source,
                    target_qualities = ConScape.sparse(hab_qual_target),
                    costs = ConScape.mapnz(x -> -log(x), adjacency_matrix))

coarse_target_qualities = ConScape.coarse_graining(g, land_mark)

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

end
