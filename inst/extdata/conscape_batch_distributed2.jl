using Logging
using ConScape
using SparseArrays
using Statistics

@everywhere using Logging, ConScape, SparseArrays, Statistics, SharedArrays

function conscape_batch_distributed(src_dir, mov_dir, target_dir, out_dir,
                                    r_targets::Vector{String}, r_sources::Vector{String},
                                    r_res::Vector{String}, land_mark, theta, exp_d, NA_val,
                                    max_retries::Int, progress::Bool)

    n = length(r_targets)
    results = Vector{Union{String, Nothing}}(undef, n)

    progress_counter = SharedArray{Int}(1)
    progress_counter[1] = 0
    print_lock = ReentrantLock()

    inputs = collect(1:n)  # just task indices

    results = pmap(inputs) do i
        iter = "-iter_$i"
        attempt = 1
        last_error = ""

        while attempt ≤ max_retries
            try
                result = conscape(src_dir, mov_dir, target_dir, out_dir,
                                  r_targets[i], r_sources[i], r_res[i],
                                  land_mark, theta, exp_d, NA_val, iter)

                if result == "done"
                    if progress
                        lock(print_lock) do
                            progress_counter[1] += 1
                            pct = round(Int, 100 * progress_counter[1] / n)
                            println(" $pct% completed")
                        end
                    end
                    return "ok"
                else
                    last_error = "conscape returned error"
                end

            catch e
                last_error = sprint((io) -> showerror(io, e, catch_backtrace()))
            end
            attempt += 1
        end

        lock(print_lock) do
            @warn "Failed iteration $i after $max_retries attempt(s): $last_error"
        end

        return last_error
    end

    return results
end
