using Distributed
using Logging
using ConScape
using SparseArrays
using Statistics

@everywhere begin
    using Logging
    using ConScape
    using SparseArrays
    using Statistics
    Logging.global_logger(Logging.NullLogger())
end

# NOTE: max_retries is kept for backward compatibility but no retries are performed.
@everywhere function safe_conscape(i, src_dir, mov_dir, target_dir, out_dir,
                                   r_targets, r_sources, r_res,
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
                                   sensitivity_require_landmark_one)
    iter = "-iter_$i"
    try
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                conscape(src_dir, mov_dir, target_dir, out_dir,
                         r_targets[i], r_sources[i], r_res[i],
                         land_mark, theta, exp_d, NA_val, iter,
                         metrics, connectivity_function_name,
                         cost_function_name, sensitivity_wrt,
                         sensitivity_method,
                         sensitivity_landscape_measure,
                         sensitivity_unitless, sensitivity_one_out_of,
                         sensitivity_diagvalue,
                         sensitivity_target_equal_source,
                         sensitivity_require_landmark_one)
            end
        end
        return "ok"
    catch e
        msg = sprint((io) -> showerror(io, e, catch_backtrace()))
        println(msg)
        return msg
    end
end

function conscape_batch_distributed(src_dir, mov_dir, target_dir, out_dir,
                                    r_targets::Vector{String}, r_sources::Vector{String},
                                    r_res::Vector{String}, land_mark, theta, exp_d, NA_val,
                                    max_retries::Int = 4, progress::Bool = true,
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

    # Align filenames across target/source/movement by basename.
    # This prevents accidental mixing of tiles when some files have been deleted.
    bn = sort!(collect(intersect(intersect(Set(r_targets), Set(r_sources)), Set(r_res))))
    r_targets = bn
    r_sources = bn
    r_res     = bn

    n = length(r_targets)
    results = Vector{Union{String,Nothing}}(undef, n)
    progress_counter = Base.Threads.Atomic{Int}(0)

    @sync for i in 1:n
        @async begin
            worker = workers()[(i - 1) % nworkers() + 1]
            results[i] = remotecall_fetch(safe_conscape, worker,
                                           i, src_dir, mov_dir, target_dir, out_dir,
                                           r_targets, r_sources, r_res,
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
                                           sensitivity_require_landmark_one)
            if progress
                done = Base.Threads.atomic_add!(progress_counter, 1) + 1
                pct = round(Int, 100 * done / n)
                println("Progress: $pct% completed")
            end
        end
    end

    return results
end
