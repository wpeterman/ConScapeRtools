using Base.Threads
using Logging
using ConScape
using SparseArrays
using Statistics

function conscape_batch(src_dir, mov_dir, target_dir, out_dir,
                        r_targets::Vector{String}, r_sources::Vector{String},
                        r_res::Vector{String}, land_mark, theta, exp_d, NA_val,
                        max_retries::Int, progress::Bool)

    # results = Vector{Union{String,Nothing}}(undef, length(r_targets))
    print_lock = ReentrantLock()
    progress_counter = Ref(0)

    @threads for i in 1:length(r_targets)
        iter = "-iter_$i"
        attempt = 1
        success = false

        while attempt ≤ max_retries && !success
            try
                Logging.with_logger(Logging.NullLogger()) do
                    redirect_stdout(devnull) do
                        redirect_stderr(devnull) do
                            conscape(src_dir, mov_dir, target_dir, out_dir,
                                     r_targets[i], r_sources[i], r_res[i],
                                     land_mark, theta, exp_d, NA_val, iter)
                        end
                    end
                end

                # results[i] = "ok"
                success = true

                if progress
                    lock(print_lock) do
                        progress_counter[] += 1
                        pct = round(Int, 100 * progress_counter[] / length(r_targets))
                        println("Progress: $pct% completed")
                    end
                end

            catch e
                if attempt == max_retries
                    # results[i] = sprint((io) -> showerror(io, e, catch_backtrace()))
                    lock(print_lock) do
                        @warn "Failed iteration $i after $attempt attempt(s)" exception=(e, catch_backtrace())
                    end
                end
                attempt += 1
            end
        end
    end

    # return results
    return "Finished processing tiles"
end
