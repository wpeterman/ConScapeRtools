# SPDX-License-Identifier: GPL-2.0-only
# _paired_metric_surfaces and _streamed_default_surfaces adapt numerical routines
# by the ConScape.jl contributors (ConScape 0.3.0, GPL-2). See COPYRIGHTS in the
# installed package for source files, source versions, and retained notices.
# Modified by the ConScapeRtools authors, 2026-10-08: shared proximity, bounded
# target-block solves, sparse factorization reuse, and accumulated contributions.

using Logging
using ConScape
using LinearAlgebra
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
    mkpath(out_dir)
    for (_, (dir_name, _)) in specs
        outdir = joinpath(out_dir, dir_name)
        mkpath(outdir)
    end
end

function _write_surface(out_dir, specs, id, iter, matrix, meta)
    dir_name, prefix = specs[id]
    ConScape.writeasc(joinpath(out_dir, dir_name, prefix * iter * ".asc"), matrix, meta)
end

function _positive_count(x)
    if x isa SparseMatrixCSC
        return nnz(x)
    else
        return count(v -> isfinite(v) && v > 0.0, x)
    end
end

function _gridsrp_failure_message(e, mov_prob, hab_qual_source, coarse_target_qualities, theta)
    msg = sprint(showerror, e)
    movement_cells = _positive_count(mov_prob)
    source_cells = _positive_count(hab_qual_source)
    target_cells = _positive_count(coarse_target_qualities)
    context = "movement cells > 0: $movement_cells; source cells > 0: $source_cells; coarse target cells > 0: $target_cells; theta: $theta"

    if e isa SingularException || occursin("SingularException", msg)
        return "GridRSP failed because ConScape built a singular linear system. " *
               "This usually indicates a disconnected or degenerate movement graph, " *
               "zero or isolated source/target cells after masking and coarse graining, " *
               "or an extreme theta value for this landscape. $context. Original error: $msg"
    elseif occursin("reducing over an empty collection", msg)
        return "GridRSP failed because ConScape found no usable graph entries after masking, " *
               "cost transformation, and theta weighting. Check that affinities are positive " *
               "where movement is possible and that source/target cells overlap reachable habitat. " *
               "$context. Original error: $msg"
    else
        return "GridRSP failed: $msg. $context"
    end
end

function _with_conscape_logging(f::Function, quiet::Bool)
    # Logger scope is task-local. Stream redirection is process-wide and unsafe
    # when tiles run on different Julia threads. Numerical warnings stay visible.
    Logging.with_logger(Logging.ConsoleLogger(stderr, quiet ? Logging.Warn : Logging.Info)) do
        f()
    end
end

# Reuse the proximity matrix when both default metrics are requested. These
# calls use ConScape's own kernels, with exactly its target ordering and no
# change to weights or diagonal treatment. No extra graph-sized cache survives.
function _paired_metric_surfaces(h, conn_fun, distance_transformation)
    proximity = conn_fun(h)
    if conn_fun <: ConScape.DistanceFunction
        map!(distance_transformation, proximity, proximity)
    end
    fcon = ConScape.connected_habitat(h, proximity)
    targetidx, targetnodes = ConScape._targetidx_and_nodes(h.g)
    coordinates = h.g.id_to_grid_coordinate_list
    bet = ConScape.RSP_betweenness_kweighted(
        h.W, h.Z,
        [h.g.source_qualities[i] for i in coordinates],
        [h.g.target_qualities[i] for i in intersect(coordinates, targetidx)],
        proximity, targetnodes)
    btwn = fill(NaN, h.g.nrows, h.g.ncols)
    btwn[coordinates] = bet
    return fcon, btwn
end

# The two default metrics are sums over targets. Solve their columns in bounded
# blocks, reusing one sparse factorization for forward and transposed solves.
# This is the same full-graph calculation, not an additional spatial
# approximation. Only floating-point reduction order differs from ConScape.
function _streamed_default_surfaces(g, theta, distance_transformation;
                                    connected = true, between = true,
                                    blocksize = 128)
    blocksize > 0 || throw(ArgumentError("blocksize must be positive"))
    coordinates = g.id_to_grid_coordinate_list
    # Preserve the upstream node order, including explicitly stored sparse zeros.
    targetidx = Set(ConScape._targetidx(g.target_qualities, coordinates))
    targetnodes = findall(in(targetidx), coordinates)
    isempty(targetnodes) && throw(ArgumentError("No retained target nodes"))
    source = [g.source_qualities[i] for i in coordinates]
    target = [g.target_qualities[coordinates[i]] for i in targetnodes]
    W = ConScape._W(ConScape._Pref(g.affinities), theta, g.costmatrix)
    factor = lu(I - W)
    CW = g.costmatrix .* W
    n = length(coordinates)
    fcon = zeros(n)
    btwn = zeros(n)
    tiny_z = false
    min_edge = minimum(nonzeros(CW))
    for first in 1:blocksize:length(targetnodes)
        cols = first:min(first + blocksize - 1, length(targetnodes))
        nodes = targetnodes[cols]
        rhs = zeros(n, length(cols))
        for (j, node) in enumerate(nodes)
            rhs[node, j] = 1.0
        end
        Z = factor \ rhs
        tiny_z |= minimum(Z) * min_edge == 0
        # Expected cost: (I-W)^-1 (C .* W) Z, divided by Z, then
        # subtract each target's diagonal value before transforming distance.
        K = factor \ (CW * Z)
        K ./= Z
        replace!(K, NaN => Inf)
        diagonal = [K[node, j] for (j, node) in enumerate(nodes)]
        K .-= diagonal'
        map!(distance_transformation, K, K)
        qt = target[cols]
        if connected
            fcon .+= source .* (K * qt)
        end
        if between
            K .*= source .* qt'
            # Scaling is local to a block and cancels after its linear solve.
            # Keep upstream overflow protection and inverse-zero convention.
            scale = max(1.0, maximum(K))
            totals = vec(sum(K, dims = 1)) ./ scale
            for j in axes(K, 2), i in axes(K, 1)
                zi = inv(Z[i, j])
                zi = isfinite(zi) ? zi : floatmax(Float64)
                K[i, j] *= inv(scale) * zi
            end
            for (j, node) in enumerate(nodes)
                zi = inv(Z[node, j])
                K[node, j] -= totals[j] * (isfinite(zi) ? zi : floatmax(Float64))
            end
            B = factor' \ K
            B .*= scale .* Z
            btwn .+= vec(sum(B, dims = 2))
        end
    end
    if tiny_z
        @warn "Z-matrix contains too small values, which can lead to inaccurate results. Check graph connectivity or decrease theta."
    end
    function surface(values)
        out = fill(NaN, g.nrows, g.ncols)
        out[coordinates] = values
        out
    end
    return connected ? surface(fcon) : nothing, between ? surface(btwn) : nothing
end

function _check_retained_habitat(g, movement, source, target)
    retained = falses(size(movement))
    retained[g.id_to_grid_coordinate_list] .= true
    dropped = count(eachindex(movement)) do i
        !retained[i] && movement[i] > 0 && (source[i] > 0 || target[i] > 0)
    end
    dropped == 0 || throw(ArgumentError(
        "ConScape's largest-component pruning would discard $dropped positive-quality cells. " *
        "Analyze disconnected components separately or revise the movement mask explicitly."))
    nothing
end

function _check_coarse_target_mass(g, coarse, landmark)
    landmark == 1 && return nothing
    coordinates = g.id_to_grid_coordinate_list
    original_mass = sum(g.target_qualities[i] for i in coordinates)
    coarse_mass = sum(coarse[i] for i in coordinates)
    original_mass > 0 && coarse_mass == 0 && throw(ArgumentError(
        "Coarse graining removed all reachable target quality. Use landmark = 1 or revise the grid alignment."))
    if !isapprox(original_mass, coarse_mass; rtol = 1e-10, atol = 0)
        @warn "Coarse graining changed reachable target-quality total; use landmark = 1 to avoid masked-landmark and partial-block effects" landmark original_mass coarse_mass
    end
    nothing
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
                  sensitivity_require_landmark_one = true,
                  quiet = true)

    _with_conscape_logging(quiet) do
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
                            throw(ErrorException("ConScape.sensitivity is not available in the installed ConScape. Run conscape_sensitivity_setup(jl_home) from R to install the ConScape sensitivity branch, then retry."))
                        elseif sensitivity_method == "simulation" && !isdefined(ConScape, :sensitivity_simulation)
                            throw(ErrorException("ConScape.sensitivity_simulation is not available in the installed ConScape. Run conscape_sensitivity_setup(jl_home) from R to install the ConScape sensitivity branch, then retry."))
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

                    for (label, values) in (("affinities", mov_prob),
                                             ("source qualities", hab_qual_source),
                                             ("target qualities", hab_qual_target))
                        any(v -> !isnan(v) && (!isfinite(v) || v < 0), values) &&
                            throw(ArgumentError("$label must contain non-negative finite values or NA"))
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
                    _check_retained_habitat(g, mov_prob, hab_qual_source, hab_qual_target)

                    coarse_target_qualities = try
                        # npix=1 is exactly the original target matrix. Skipping
                        # coarse_graining also avoids its repeated node searches.
                        land_mark == 1 ? g.target_qualities : ConScape.coarse_graining(g, land_mark)
                    catch e
                        throw(ErrorException("coarse_graining failed: $e"))
                    end
                    _check_coarse_target_mass(g, coarse_target_qualities, land_mark)

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
                        land_mark == 1 ? g : ConScape.Grid(size(mov_prob)...,
                                      affinities = adjacency_matrix,
                                      source_qualities = hab_qual_source,
                                      target_qualities = coarse_target_qualities,
                                      costs = costs)
                    catch e
                        throw(ErrorException("coarse Grid construction failed: $e"))
                    end

                    distance_transformation = x -> exp(-x / exp_d)
                    if conn_fun === ConScape.expected_cost && isempty(sensitivity_wrt) &&
                       all(in(("connected_habitat", "betweenness_kweighted")), metric_ids)
                        try
                            fcon, btwn = _streamed_default_surfaces(g_coarse, theta,
                                distance_transformation;
                                connected = "connected_habitat" in metric_ids,
                                between = "betweenness_kweighted" in metric_ids)
                            fcon === nothing || _write_surface(out_dir, specs, "connected_habitat", iter, fcon, meta_p)
                            btwn === nothing || _write_surface(out_dir, specs, "betweenness_kweighted", iter, btwn, meta_p)
                        catch e
                            throw(ErrorException(_gridsrp_failure_message(e, mov_prob,
                                hab_qual_source, coarse_target_qualities, theta)))
                        end
                        return "done"
                    end

                    h_coarse = try
                        ConScape.GridRSP(g_coarse, θ = theta)
                    catch e
                        try
                            # Only retry a keyword-API mismatch, not a failed solve.
                            e isa MethodError || rethrow()
                            ConScape.GridRSP(g_coarse, theta = theta)
                        catch e2
                            throw(ErrorException(_gridsrp_failure_message(
                                e2,
                                mov_prob,
                                hab_qual_source,
                                coarse_target_qualities,
                                theta)))
                        end
                    end

                    distance_transformation = x -> exp(-x / exp_d)

                    try
                        paired = "connected_habitat" in metric_ids && "betweenness_kweighted" in metric_ids
                        if paired
                            fcon, btwn = _paired_metric_surfaces(h_coarse, conn_fun, distance_transformation)
                            _write_surface(out_dir, specs, "connected_habitat", iter, fcon, meta_p)
                            _write_surface(out_dir, specs, "betweenness_kweighted", iter, btwn, meta_p)
                        end
                        for metric in metric_ids
                            if paired && metric in ("connected_habitat", "betweenness_kweighted")
                                continue
                            end
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

    return "done"
end
