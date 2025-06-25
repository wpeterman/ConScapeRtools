using Logging

function conscape(src_dir, mov_dir, target_dir, out_dir, r_target, r_source, r_res,
                  land_mark, theta, exp_d, NA_val, iter)

    Logging.with_logger(Logging.NullLogger()) do
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                try
                    targetdir = joinpath(target_dir)
                    movdir = joinpath(mov_dir)
                    srcdir = joinpath(src_dir)
                    outdir_btwn = joinpath(out_dir, "btwn")
                    outdir_fcon = joinpath(out_dir, "fcon")

                    if !isdir(out_dir); mkdir(out_dir); end
                    if !isdir(outdir_btwn); mkdir(outdir_btwn); end
                    if !isdir(outdir_fcon); mkdir(outdir_fcon); end

                    hab_qual_target, meta_t = try
                        ConScape.readasc(joinpath(targetdir, r_target))
                    catch e
                        throw(ErrorException("readasc r_target failed: $e"))
                    end

                    hab_qual_source, meta_s = try
                        ConScape.readasc(joinpath(srcdir, r_source))
                    catch e
                        throw(ErrorException("readasc r_source failed: $e"))
                    end

                    mov_prob, meta_p = try
                        ConScape.readasc(joinpath(movdir, r_res))
                    catch e
                        throw(ErrorException("readasc r_res failed: $e"))
                    end

                    try
                        non_matches = findall(xor.(isnan.(mov_prob), isnan.(hab_qual_target)))
                        mov_prob[non_matches] .= 0
                        hab_qual_target[non_matches] .= NA_val

                        non_matches = findall(xor.(isnan.(mov_prob), isnan.(hab_qual_source)))
                        mov_prob[non_matches] .= NA_val
                        hab_qual_source[non_matches] .= NA_val
                    catch e
                        throw(ErrorException("NA handling failed: $e"))
                    end

                    adjacency_matrix = try
                        ConScape.graph_matrix_from_raster(mov_prob)
                    catch e
                        throw(ErrorException("graph_matrix_from_raster failed: $e"))
                    end

                    g = try
                        ConScape.Grid(size(mov_prob)...,
                                      affinities = adjacency_matrix,
                                      source_qualities = hab_qual_source,
                                      target_qualities = ConScape.sparse(hab_qual_target),
                                      costs = ConScape.mapnz(x -> -log(x), adjacency_matrix))
                    catch e
                        throw(ErrorException("Grid construction failed: $e"))
                    end

                    coarse_target_qualities = try
                        ConScape.coarse_graining(g, land_mark)
                    catch e
                        throw(ErrorException("coarse_graining failed: $e"))
                    end

                    g_coarse = try
                        ConScape.Grid(size(mov_prob)...,
                                      affinities = adjacency_matrix,
                                      source_qualities = hab_qual_source,
                                      target_qualities = coarse_target_qualities,
                                      costs = ConScape.mapnz(x -> -log(x), adjacency_matrix))
                    catch e
                        throw(ErrorException("coarse Grid construction failed: $e"))
                    end

                    h_coarse = try
                        ConScape.GridRSP(g_coarse, θ = theta)
                    catch e
                        throw(ErrorException("GridRSP failed: $e"))
                    end

                    func_con = try
                        ConScape.connected_habitat(h_coarse,
                            distance_transformation = x -> exp(-x / exp_d))
                    catch e
                        throw(ErrorException("connected_habitat failed: $e"))
                    end

                    kbetw = try
                        ConScape.betweenness_kweighted(h_coarse,
                            distance_transformation = x -> exp(-x / exp_d))
                    catch e
                        throw(ErrorException("betweenness_kweighted failed: $e"))
                    end

                    try
                        ConScape.writeasc(joinpath(outdir_btwn, "betweenness" * iter * ".asc"), kbetw, meta_p)
                        ConScape.writeasc(joinpath(outdir_fcon, "fcon" * iter * ".asc"), func_con, meta_p)
                    catch e
                        throw(ErrorException("writeasc failed: $e"))
                    end

                catch e
                    throw(ErrorException("Unexpected failure in conscape: $e"))
                end
            end
        end
    end

    return "done"
end
