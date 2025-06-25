using Base.Threads
using Logging
using ConScape
using SparseArrays
using Statistics

function conscape_batch(src_dir, mov_dir, target_dir, out_dir,
                        r_targets::Vector{String}, r_sources::Vector{String},
                        r_res::Vector{String}, land_mark, theta, exp_d, NA_val,
                        max_retries::Int, progress::Bool)

    results = Vector{Union{String, Nothing}}(undef, length(r_targets))
    print_lock = ReentrantLock()
    progress_counter = Ref(0)

    @threads for i in 1:length(r_targets)
        iter = "-iter_$i"
        attempt = 1
        success = false
        last_error = ""

        while attempt ≤ max_retries && !success
            try
                result = conscape(src_dir, mov_dir, target_dir, out_dir,
                                  r_targets[i], r_sources[i], r_res[i],
                                  land_mark, theta, exp_d, NA_val, iter)

                if result == "done"
                    success = true
                    results[i] = "ok"

                    if progress
                        lock(print_lock) do
                            progress_counter[] += 1
                            pct = round(Int, 100 * progress_counter[] / length(r_targets))
                            println("Progress: $pct% completed")
                        end
                    end
                else
                    last_error = "conscape returned error"
                end

            catch e
                last_error = sprint((io) -> showerror(io, e, catch_backtrace()))
            end
            attempt += 1
        end

        if !success
            lock(print_lock) do
                @warn "Failed iteration $i after $max_retries attempt(s): $last_error"
            end
            results[i] = last_error
        end
    end

    return results
end
