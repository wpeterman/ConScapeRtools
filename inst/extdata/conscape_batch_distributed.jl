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
                                   max_retries)
    iter = "-iter_$i"
    try
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                conscape(src_dir, mov_dir, target_dir, out_dir,
                         r_targets[i], r_sources[i], r_res[i],
                         land_mark, theta, exp_d, NA_val, iter)
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
                                    max_retries::Int = 4, progress::Bool = true)

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
                                          max_retries)
            if progress
                done = Base.Threads.atomic_add!(progress_counter, 1) + 1
                pct = round(Int, 100 * done / n)
                println("Progress: $pct% completed")
            end
        end
    end

    return results
end
