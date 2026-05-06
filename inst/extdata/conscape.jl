using Logging
using ConScape
using SparseArrays
using Statistics

function _as_string_vector(x)
    if x === nothing
        return String[]
    elseif x isa AbstractVector
        return filter(!isempty, String.(x))
    elseif x isa AbstractString
        value = String(x)
        return isempty(value) ? String[] : [value]
    else
        return String[]
    end
end

function _costs_for_grid(cost_function, adjacency_matrix)
    cost_function = lowercase(strip(String(cost_function)))
    if cost_function in ["minuslog", "-log", "x -> -log(x)"]
        return ConScape.MinusLog()
    elseif cost_function in ["inverse", "inv", "x -> inv(x)"]
        return ConScape.Inv()
    elseif cost_function == "odds_against"
        return ConScape.OddsAgainst()
    elseif cost_function == "odds_for"
        return ConScape.OddsFor()
    elseif cost_function == "expminus"
        return ConScape.ExpMinus()
    elseif occursin("->", cost_function)
        return ConScape.mapnz(eval(Meta.parse(cost_function)), adjacency_matrix)
    else
        throw(ArgumentError("unsupported cost_function: $cost_function"))
    end
end

function _connectivity_function(name)
    name = lowercase(strip(String(name)))
    if name == "expected_cost"
        return ConScape.expected_cost
    elseif name == "least_cost_distance"
        return ConScape.least_cost_distance
    elseif name == "free_energy_distance"
        return ConScape.free_energy_distance
    elseif name == "survival_probability"
        return ConScape.survival_probability
    elseif name == "power_mean_proximity"
        return ConScape.power_mean_proximity
    else
        throw(ArgumentError("unsupported connectivity_function: $name"))
    end
end

function _metric_specs(metrics, sensitivity_wrt, sensitivity_unitless)
    specs = Dict{String,Tuple{String,String}}()

    for metric in _as_string_vector(metrics)
        if metric == "connected_habitat"
            specs["connected_habitat"] = ("fcon", "fcon")
        elseif metric == "betweenness_kweighted"
            specs["betweenness_kweighted"] = ("btwn", "betweenness")
        elseif metric == "betweenness_qweighted"
            specs["betweenness_qweighted"] = ("btwn_qweighted", "betweenness_qweighted")
        elseif metric == "criticality"
            specs["criticality"] = ("criticality", "criticality")
        else
            throw(ArgumentError("unsupported metric: $metric"))
        end
    end

    prefix = sensitivity_unitless ? "elasticity" : "sensitivity"
    for wrt in _as_string_vector(sensitivity_wrt)
        label = if wrt == "Q"
            "quality"
        elseif wrt == "A"
            "affinity"
        elseif wrt == "C"
            "cost"
        elseif wrt == "A&C=f(A)"
            "affinity_cost_linked"
        elseif wrt == "C&A=f(C)"
            "cost_affinity_linked"
        else
            throw(ArgumentError("unsupported sensitivity wrt: $wrt"))
        end
        id = prefix * "_" * label
        specs[id] = (id, id)
    end

    return specs
end

function _ensure_output_dirs(out_dir, specs)
    if !isdir(out_dir); mkdir(out_dir); end
    for (_, (dir_name, _)) in specs
        outdir = joinpath(out_dir, dir_name)
        if !isdir(outdir); mkdir(outdir); end
    end
end

function _write_surface(out_dir, specs, id, iter, matrix, meta)
    dir_name, prefix = specs[id]
    ConScape.writeasc(joinpath(out_dir, dir_name, prefix * iter * ".asc"), matrix, meta)
end

function conscape(src_dir, mov_dir, target_dir, out_dir, r_target, r_source, r_res,
                  land_mark, theta, exp_d, NA_val, iter,
                  metrics = ["betweenness_kweighted", "connected_habitat"],
                  connectivity_function_name = "expected_cost",
                  cost_function_name = "minuslog",
                  sensitivity_wrt = String[],
                  sensitivity_method = "analytical",
                  sensitivity_landscape_measure = "sum",
                  sensitivity_unitless = true,
                  sensitivity_one_out_of = 1,
                  sensitivity_diagvalue = nothing,
                  sensitivity_target_equal_source = true,
                  sensitivity_require_landmark_one = true)

    Logging.with_logger(Logging.NullLogger()) do
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                try
                    targetdir = joinpath(target_dir)
                    movdir = joinpath(mov_dir)
                    srcdir = joinpath(src_dir)

                    sensitivity_wrt = _as_string_vector(sensitivity_wrt)
                    metric_ids = _as_string_vector(metrics)
                    specs = _metric_specs(metric_ids, sensitivity_wrt, sensitivity_unitless)
                    _ensure_output_dirs(out_dir, specs)

                    if !isempty(sensitivity_wrt)
                        if sensitivity_require_landmark_one && land_mark != 1
                            throw(ErrorException("Sensitivity requires land_mark = 1 so target qualities remain equal to source qualities."))
                        end
                        if sensitivity_method == "analytical" && !isdefined(ConScape, :sensitivity)
                            throw(ErrorException("ConScape.sensitivity is not available. Install the ConScape sensitivity branch or a release that includes it."))
                        elseif sensitivity_method == "simulation" && !isdefined(ConScape, :sensitivity_simulation)
                            throw(ErrorException("ConScape.sensitivity_simulation is not available. Install the ConScape sensitivity branch or a release that includes it."))
                        end
                    end

                    conn_fun = _connectivity_function(connectivity_function_name)

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
                        mov_nan = isnan.(mov_prob)

                        if any(mov_nan)
                            mov_prob[mov_nan] .= 0.0
                            hab_qual_source[mov_nan] .= 0.0
                            hab_qual_target[mov_nan] .= 0.0
                        end

                        hab_qual_target[isnan.(hab_qual_target)] .= 0.0
                        hab_qual_source[isnan.(hab_qual_source)] .= 0.0
                    catch e
                        throw(ErrorException("NA handling failed: $e"))
                    end

                    adjacency_matrix = try
                        ConScape.graph_matrix_from_raster(mov_prob)
                    catch e
                        throw(ErrorException("graph_matrix_from_raster failed: $e"))
                    end

                    costs = _costs_for_grid(cost_function_name, adjacency_matrix)

                    g = try
                        ConScape.Grid(size(mov_prob)...,
                                      affinities = adjacency_matrix,
                                      source_qualities = hab_qual_source,
                                      target_qualities = ConScape.sparse(hab_qual_target),
                                      costs = costs)
                    catch e
                        throw(ErrorException("Grid construction failed: $e"))
                    end

                    coarse_target_qualities = try
                        ConScape.coarse_graining(g, land_mark)
                    catch e
                        throw(ErrorException("coarse_graining failed: $e"))
                    end

                    has_targets(x) = x isa SparseMatrixCSC ? (nnz(x) > 0) : any(x .> 0.0)

                    if !has_targets(coarse_target_qualities)
                        blank = fill(NaN, size(mov_prob)...)
                        try
                            for id in keys(specs)
                                _write_surface(out_dir, specs, id, iter, blank, meta_p)
                            end
                        catch e
                            throw(ErrorException("writeasc (blank outputs) failed: $e"))
                        end
                        return "skipped_empty_targets"
                    end

                    g_coarse = try
                        ConScape.Grid(size(mov_prob)...,
                                      affinities = adjacency_matrix,
                                      source_qualities = hab_qual_source,
                                      target_qualities = coarse_target_qualities,
                                      costs = costs)
                    catch e
                        throw(ErrorException("coarse Grid construction failed: $e"))
                    end

                    h_coarse = try
                        ConScape.GridRSP(g_coarse, theta = theta)
                    catch e
                        try
                            ConScape.GridRSP(g_coarse, θ = theta)
                        catch e2
                            throw(ErrorException("GridRSP failed: $e2"))
                        end
                    end

                    distance_transformation = x -> exp(-x / exp_d)

                    try
                        for metric in metric_ids
                            if metric == "connected_habitat"
                                value = ConScape.connected_habitat(
                                    h_coarse,
                                    connectivity_function = conn_fun,
                                    distance_transformation = distance_transformation)
                                _write_surface(out_dir, specs, metric, iter, value, meta_p)
                            elseif metric == "betweenness_kweighted"
                                value = ConScape.betweenness_kweighted(
                                    h_coarse,
                                    connectivity_function = conn_fun,
                                    distance_transformation = distance_transformation)
                                _write_surface(out_dir, specs, metric, iter, value, meta_p)
                            elseif metric == "betweenness_qweighted"
                                value = ConScape.betweenness_qweighted(h_coarse)
                                _write_surface(out_dir, specs, metric, iter, value, meta_p)
                            elseif metric == "criticality"
                                value = ConScape.criticality(
                                    h_coarse,
                                    distance_transformation = distance_transformation)
                                _write_surface(out_dir, specs, metric, iter, value, meta_p)
                            end
                        end
                    catch e
                        throw(ErrorException("metric computation failed: $e"))
                    end

                    if !isempty(sensitivity_wrt)
                        alpha = 1 / exp_d
                        try
                            for wrt in sensitivity_wrt
                                id = (sensitivity_unitless ? "elasticity_" : "sensitivity_") *
                                     (wrt == "Q" ? "quality" :
                                      wrt == "A" ? "affinity" :
                                      wrt == "C" ? "cost" :
                                      wrt == "A&C=f(A)" ? "affinity_cost_linked" :
                                      "cost_affinity_linked")

                                value = if sensitivity_method == "simulation"
                                    ConScape.sensitivity_simulation(
                                        h_coarse,
                                        connectivity_function = conn_fun,
                                        distance_transformation = ConScape.ExpMinus(),
                                        α = alpha,
                                        wrt = wrt,
                                        landscape_measure = sensitivity_landscape_measure,
                                        unitless = sensitivity_unitless,
                                        diagvalue = sensitivity_diagvalue,
                                        target_equal_source = sensitivity_target_equal_source,
                                        one_out_of = sensitivity_one_out_of)
                                else
                                    ConScape.sensitivity(
                                        h_coarse,
                                        connectivity_function = conn_fun,
                                        distance_transformation = ConScape.ExpMinus(),
                                        α = alpha,
                                        wrt = wrt,
                                        landscape_measure = sensitivity_landscape_measure,
                                        unitless = sensitivity_unitless,
                                        diagvalue = sensitivity_diagvalue,
                                        target_equal_source = sensitivity_target_equal_source)
                                end
                                _write_surface(out_dir, specs, id, iter, value, meta_p)
                            end
                        catch e
                            throw(ErrorException("sensitivity computation failed: $e"))
                        end
                    end

                catch e
                    throw(ErrorException("Unexpected failure in conscape: $e"))
                end
            end
        end
    end

    return "done"
end
