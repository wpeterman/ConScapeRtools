using Distributed
using LinearAlgebra
using Logging
using ConScape
using SparseArrays
using Statistics

# The caller includes this file only on its own worker pool. Do not broadcast
# definitions or imports to unrelated workers in a reused Julia session.
function _set_blas_threads(n::Integer)
    try
        BLAS.set_num_threads(max(1, n))
    catch e
        @warn "Could not set BLAS threads" exception=(e, catch_backtrace())
    end
    return nothing
end

function _csv_escape(x)
    s = replace(String(x), "\"" => "\"\"")
    return "\"" * s * "\""
end

function _write_batch_diagnostics(out_dir, rows)
    path = joinpath(out_dir, "conscape_batch_diagnostics.csv")
    open(path, "w") do io
        println(io, "index,tile,status,attempts,seconds,work_estimate,error")
        for row in rows
            println(io, join((
                row.index,
                _csv_escape(row.tile),
                _csv_escape(row.status),
                row.attempts,
                row.seconds,
                row.work_estimate,
                _csv_escape(row.error)
            ), ","))
        end
    end
    return path
end

# NOTE: max_retries is kept for backward compatibility but no retries are performed.
function safe_conscape(i, src_dir, mov_dir, target_dir, out_dir,
                                   r_target, r_source, r_move,
                                   land_mark, theta, exp_d, NA_val,
                                   max_retries,
                                   metrics,
                                   connectivity_function_name,
                                   cost_function_name,
                                   sensitivity_wrt,
                                   sensitivity_method,
                                   sensitivity_landscape_measure,
                                   sensitivity_unitless,
                                   sensitivity_one_out_of,
                                   sensitivity_diagvalue,
                                   sensitivity_target_equal_source,
                                   sensitivity_require_landmark_one,
                                   blas_threads)
    old_blas = BLAS.get_num_threads()
    _set_blas_threads(blas_threads)
    iter = "-" * splitext(r_target)[1]
    started = time()
    try
        result = conscape(src_dir, mov_dir, target_dir, out_dir,
                          r_target, r_source, r_move,
                          land_mark, theta, exp_d, NA_val, iter,
                          metrics, connectivity_function_name,
                          cost_function_name, sensitivity_wrt,
                          sensitivity_method,
                          sensitivity_landscape_measure,
                          sensitivity_unitless, sensitivity_one_out_of,
                          sensitivity_diagvalue,
                          sensitivity_target_equal_source,
                          sensitivity_require_landmark_one,
                          true)
        status = result in ("done", "skipped_empty_targets") ? "ok" : "failed"
        return (
            index = i,
            tile = r_target,
            status = status,
            attempts = 1,
            seconds = round(time() - started; digits = 3),
            error = status == "ok" ? "" : "conscape returned $result"
        )
    catch e
        msg = sprint((io) -> showerror(io, e, catch_backtrace()))
        println(msg)
        return (
            index = i,
            tile = r_target,
            status = "failed",
            attempts = 1,
            seconds = round(time() - started; digits = 3),
            error = msg
        )
    finally
        BLAS.set_num_threads(old_blas)
    end
end

function conscape_batch_distributed(src_dir, mov_dir, target_dir, out_dir,
                                    r_targets, r_sources,
                                    r_res, land_mark, theta, exp_d, NA_val,
                                    max_retries::Int = 4, progress::Bool = true,
                                    blas_threads::Int = 1,
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
                                    worker_ids = workers())

    worker_ids = worker_ids isa Integer ? [Int(worker_ids)] : Int.(worker_ids)
    isempty(worker_ids) && throw(ArgumentError("No distributed workers supplied"))

    # r_targets/r_sources/r_res arrive via JuliaConnectoR from R character
    # vectors. R does not distinguish a length-1 vector from a scalar, so
    # JuliaConnectoR auto-unboxes a single-tile design's filename list into a
    # bare Julia String instead of a Vector{String}. A strict
    # ::Vector{String} signature then rejects that call with a MethodError
    # for any landscape small enough to produce exactly one tile. Normalize
    # explicitly with _as_string_vector (defined in conscape.jl, included
    # via @everywhere before this file) so both scalar and vector inputs
    # work uniformly.
    r_targets = _as_string_vector(r_targets)
    r_sources = _as_string_vector(r_sources)
    r_res     = _as_string_vector(r_res)

    # Never silently omit an unmatched tile.
    Set(r_targets) == Set(r_sources) == Set(r_res) ||
        throw(ArgumentError("Target, source, and movement tile filenames must match"))
    bn = sort!(unique(r_targets))
    r_targets = bn
    r_sources = bn
    r_res     = bn

    n = length(r_targets)
    results = Vector{Union{String,Nothing}}(undef, n)
    diagnostics = Vector{NamedTuple}(undef, n)
    progress_counter = Base.Threads.Atomic{Int}(0)
    nextjob = Base.Threads.Atomic{Int}(1)
    work_estimates = [
        filesize(joinpath(target_dir, r_targets[i])) +
        filesize(joinpath(src_dir, r_sources[i])) +
        filesize(joinpath(mov_dir, r_res[i]))
        for i in 1:n
    ]
    job_order = sortperm(work_estimates, rev = true)

    # One outstanding call per worker bounds live graphs and transfers only
    # three scalar filenames per job instead of three complete manifests.
    @sync for worker in worker_ids
        @async begin
          while true
            job = Base.Threads.atomic_add!(nextjob, 1)
            job > n && break
            i = job_order[job]
            row = remotecall_fetch(safe_conscape, worker,
                                   i, src_dir, mov_dir, target_dir, out_dir,
                                   r_targets[i], r_sources[i], r_res[i],
                                   land_mark, theta, exp_d, NA_val,
                                   max_retries, metrics,
                                   connectivity_function_name,
                                   cost_function_name,
                                   sensitivity_wrt,
                                   sensitivity_method,
                                   sensitivity_landscape_measure,
                                   sensitivity_unitless,
                                   sensitivity_one_out_of,
                                   sensitivity_diagvalue,
                                   sensitivity_target_equal_source,
                                   sensitivity_require_landmark_one,
                                   blas_threads)
            diagnostics[i] = merge(row, (work_estimate = work_estimates[i],))
            results[i] = row.status == "ok" ? "ok" : row.error
            if progress
                done = Base.Threads.atomic_add!(progress_counter, 1) + 1
                pct = round(Int, 100 * done / n)
                println("Progress: $pct% completed")
            end
          end
        end
    end

    _write_batch_diagnostics(out_dir, diagnostics)
    return results
end
