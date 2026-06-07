using ConScape
using LinearAlgebra
using ArchGDAL
using Rasters

function _dev_require_api(mode::AbstractString)
    required = Symbol[:Problem, :WindowedProblem, :solve]
    if mode == "batch"
        append!(required, [:BatchProblem, :assess])
    end
    missing = [String(s) for s in required if !isdefined(ConScape, s)]
    if !isempty(missing)
        error("The experimental ConScape dev backend requires a ConScape installation with " *
              join(missing, ", ") * ". Install the ConScape.jl alg_efficiency/dev branch to use this backend.")
    end
    return nothing
end

function _dev_install_solver_compat_patch()
    @eval ConScape begin
        function _init_dense!(
            ws::NamedTuple,
            solver::Solver,
            cm::FundamentalMeasure,
            p::AbstractProblem,
            grid::Grid;
            verbose=false,
            reuse_output=false,
        )
            verbose && println("Retrieving measures...")
            gms = graph_measures(p)
            cf = connectivity_function(p)
            verbose && println("Defining sparse arrays...")
            verbose && println("Allocating workspaces...")
            sze = _workspace_size(solver, grid)
            Z = if hastrait(needs_inv, gms)
                if haskey(ws, :Z)
                    _reshape(ws.Z, sze)
                else
                    Matrix{eltype(grid.affinities)}(undef, sze)
                end
            else
                nothing
            end
            Zⁱ = if hastrait(needs_inv, gms)
                haskey(ws, :Zⁱ) ? _reshape(ws.Zⁱ, sze) : similar(Z)
            else
                nothing
            end
            n_workspaces = count_workspaces(p)
            n_permuted_workspaces = count_permuted_workspaces(p)
            workspaces = if haskey(ws, :workspaces)
                [_reshape(w, size(Z)) for w in ws.workspaces]
            else
                [similar(Z) for _ in 1:n_workspaces]
            end
            permuted_workspaces = if haskey(ws, :permuted_workspaces)
                [_reshape(pw, size(Z')) for pw in ws.permuted_workspaces]
            else
                [similar(Z') for _ in 1:n_permuted_workspaces]
            end
            expected_costs = if hastrait(needs_expected_cost, gms) || cf == ConScape.expected_cost
                haskey(ws, :expected_costs) ? _reshape(ws.expected_costs, size(Z)) : similar(Z)
            else
                nothing
            end
            free_energy_distances = if hastrait(needs_free_energy_distance, gms) || cf == ConScape.free_energy_distance
                haskey(ws, :free_energy_distances) ? _reshape(ws.free_energy_distances, size(Z)) : similar(Z)
            else
                nothing
            end
            proximities = if hastrait(needs_proximity, gms)
                haskey(ws, :proximities) ? _reshape(ws.proximities, size(Z)) : similar(Z)
            else
                nothing
            end
            function matrix_or_nothing(gm)
                if returntype(gm) isa ReturnsDenseSpatial
                    A = fill(NaN, size(grid))
                    A[grid.id_to_grid_coordinate_list] .= 0.0
                    A
                else
                    nothing
                end
            end
            outputs = if reuse_output && haskey(ws, :outputs)
                ws.outputs
            else
                if distance_transformation(cm) isa NamedTuple
                    map(gms) do gm
                        if needs_connectivity(gm)
                            map(distance_transformation(cm)) do dt
                                matrix_or_nothing(gm)
                            end
                        else
                            matrix_or_nothing(gm)
                        end
                    end
                else
                    map(gms) do gm
                        matrix_or_nothing(gm)
                    end
                end
            end

            subgrids = split_subgraphs(grid)
            verbose && println("Finished allocating...")
            return (; Z, Zⁱ, workspaces, permuted_workspaces,
                    free_energy_distances, expected_costs, proximities,
                    outputs, grid, g=grid, subgrids)
        end
    end
    return nothing
end
_dev_install_solver_compat_patch()

function _dev_set_blas_threads(n::Integer)
    try
        BLAS.set_num_threads(max(1, n))
    catch e
        @warn "Could not set BLAS threads" exception=(e, catch_backtrace())
    end
    return nothing
end

function _dev_cost_function(name)
    name = lowercase(strip(String(name)))
    if name == "minuslog"
        return ConScape.MinusLog()
    elseif name == "inverse"
        return ConScape.Inv()
    elseif name == "odds_against"
        return ConScape.OddsAgainst()
    elseif name == "odds_for"
        return ConScape.OddsFor()
    elseif name == "expminus"
        return ConScape.ExpMinus()
    else
        error("Unsupported cost_function for ConScape dev backend: $name")
    end
end

function _dev_measures(metrics)
    pairs = Pair{Symbol,Any}[]
    for metric in String.(metrics)
        if metric == "connected_habitat"
            push!(pairs, :fcon => ConScape.ConnectedHabitat())
        elseif metric == "betweenness_kweighted"
            push!(pairs, :btwn => ConScape.BetweennessKweighted())
        elseif metric == "betweenness_qweighted"
            push!(pairs, :btwn_qweighted => ConScape.BetweennessQweighted())
        else
            error("Unsupported metric for ConScape dev backend: $metric")
        end
    end
    return (; pairs...)
end

function _dev_problem(metrics, theta::Real, exp_d::Real, cost_function_name::AbstractString)
    graph_measures = _dev_measures(metrics)
    connectivity_measure = ConScape.ExpectedCost(
        θ = theta,
        distance_transformation = x -> exp(-x / exp_d)
    )
    return ConScape.Problem(;
        graph_measures,
        connectivity_measure,
        costs = _dev_cost_function(cost_function_name),
        solver = ConScape.VectorSolver()
    )
end

function _dev_raster_stack(target_file::String,
                           source_file::String,
                           affinities_file::String)
    paths = (
        affinities = affinities_file,
        qualities = source_file,
        target_qualities = target_file,
    )
    return Rasters.RasterStack(paths; missingval = NaN)
end

function _dev_write_outputs(result, out_dir::String, mode::String)
    output_root = joinpath(out_dir, "conscape_dev_" * mode)
    mkpath(output_root)
    written = Rasters.write(joinpath(output_root, ""), result; ext = ".tif", force = true)
    return string(written)
end

function _dev_batch_stack(dir::String, layer_names)
    pairs = Pair{Symbol,String}[]
    for name in layer_names
        path = joinpath(dir, string(name) * ".tif")
        isfile(path) || error("Missing ConScape dev batch output: $path")
        push!(pairs, Symbol(name) => path)
    end
    return Rasters.RasterStack((; pairs...); missingval = NaN)
end

function _dev_mosaic_batch(batch_problem, rast, layer_names; verbose::Bool = false)
    paths = ConScape.batch_paths(batch_problem, rast)
    dirs = filter(isdir, paths)
    isempty(dirs) && error("No completed ConScape dev batch output directories found.")
    stacks = [_dev_batch_stack(dir, layer_names) for dir in dirs]
    return Rasters.mosaic(sum, stacks; to = rast, missingval = 0.0, verbose = verbose)
end

function conscape_dev_run(target_file::String,
                          source_file::String,
                          affinities_file::String,
                          out_dir::String,
                          land_mark::Integer,
                          theta::Real,
                          exp_d::Real,
                          centersize::Integer,
                          buffer::Integer,
                          window_shape::String,
                          mode::String,
                          batch_grain,
                          batch_ext::String,
                          threaded::Bool,
                          blas_threads::Integer,
                          metrics = ["betweenness_kweighted", "connected_habitat"],
                          cost_function_name = "minuslog",
                          progress::Bool = true)
    mode = lowercase(strip(mode))
    if !(mode in ("windowed", "batch"))
        error("Unsupported ConScape dev mode: $mode")
    end
    if window_shape != "square"
        error("The current ConScape dev WindowedProblem API supports square windows; got window_shape = $window_shape")
    end
    _dev_require_api(mode)
    _dev_install_solver_compat_patch()
    _dev_set_blas_threads(blas_threads)
    if !isdir(out_dir); mkpath(out_dir); end

    rast = _dev_raster_stack(target_file, source_file, affinities_file)
    if land_mark > 1
        progress && println("Applying ConScape coarse_graining with landmark = ", land_mark)
        rast = ConScape.coarse_graining(rast, Int(land_mark))
    end
    problem = _dev_problem(metrics, theta, exp_d, cost_function_name)

    result = if mode == "windowed"
        windowed_problem = ConScape.WindowedProblem(
            problem;
            centersize = Int(centersize),
            buffer = Int(buffer),
            threaded = threaded
        )
        ConScape.solve(windowed_problem, rast; mosaic_return = true, verbose = progress)
    else
        grain_value = batch_grain === nothing ? nothing : Int(batch_grain)
        batch_problem = ConScape.BatchProblem(
            problem;
            centersize = Int(centersize),
            buffer = Int(buffer),
            datapath = joinpath(out_dir, "conscape_dev_batch_work"),
            grain = grain_value,
            ext = batch_ext
        )
        assessment = ConScape.assess(batch_problem, rast)
        progress && println("ConScape dev BatchProblem jobs: ", assessment.njobs)
        ConScape.solve(batch_problem, rast, assessment; verbose = progress)
        _dev_mosaic_batch(batch_problem, rast, keys(_dev_measures(metrics)); verbose = progress)
    end

    writer = Base.invokelatest(getfield, Main, :_dev_write_outputs)
    written = Base.invokelatest(writer, result, out_dir, mode)
    progress && println("ConScape dev ", mode, " backend wrote outputs to ", written)
    return written
end

function conscape_dev_windowed(target_file::String,
                               source_file::String,
                               affinities_file::String,
                               out_dir::String,
                               land_mark::Integer,
                               theta::Real,
                               exp_d::Real,
                               centersize::Integer,
                               buffer::Integer,
                               window_shape::String,
                               threaded::Bool,
                               blas_threads::Integer,
                               metrics = ["betweenness_kweighted", "connected_habitat"],
                               cost_function_name = "minuslog",
                               progress::Bool = true)
    runner = Base.invokelatest(getfield, Main, :conscape_dev_run)
    return Base.invokelatest(
        runner,
        target_file,
        source_file,
        affinities_file,
        out_dir,
        land_mark,
        theta,
        exp_d,
        centersize,
        buffer,
        window_shape,
        "windowed",
        nothing,
        ".tif",
        threaded,
        blas_threads,
        metrics,
        cost_function_name,
        progress
    )
end
