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

@everywhere function safe_conscape(i, src_dir, mov_dir, target_dir, out_dir,
                                   r_targets, r_sources, r_res,
                                   land_mark, theta, exp_d, NA_val,
                                   max_retries)
    iter = "-iter_$i"
    attempt = 1
    while attempt ≤ max_retries
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
            if attempt == max_retries
                return sprint((io) -> showerror(io, e, catch_backtrace()))
            end
            attempt += 1
        end
    end
end

function conscape_batch_distributed(src_dir, mov_dir, target_dir, out_dir,
                                    r_targets::Vector{String}, r_sources::Vector{String},
                                    r_res::Vector{String}, land_mark, theta, exp_d, NA_val;
                                    max_retries::Int = 4)

    n = length(r_targets)
    results = Vector{Union{String,Nothing}}(undef, n)

    @sync for i in 1:n
        @async begin
            worker = workers()[(i - 1) % nworkers() + 1]
            results[i] = remotecall_fetch(safe_conscape, worker,
                                          i, src_dir, mov_dir, target_dir, out_dir,
                                          r_targets, r_sources, r_res,
                                          land_mark, theta, exp_d, NA_val,
                                          max_retries)
        end
    end

    return results
end
