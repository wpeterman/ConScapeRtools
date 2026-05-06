using Base.Threads
using Logging
using ConScape
using SparseArrays
using Statistics

function conscape_batch(src_dir, mov_dir, target_dir, out_dir,
                        r_targets::Vector{String}, r_sources::Vector{String},
                        r_res::Vector{String}, land_mark, theta, exp_d, NA_val,
                        max_retries::Int, progress::Bool,
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
    results = Vector{Union{String, Nothing}}(undef, n)
    print_lock = ReentrantLock()
    progress_counter = Ref(0)

    @threads for i in 1:n
        # Use the filename stem (e.g., r_191) for a stable, interpretable iteration tag.
        iter = "-" * splitext(r_targets[i])[1]
        attempt = 1
        success = false
        last_error = ""

        while attempt ≤ max_retries && !success
            try
                result = conscape(src_dir, mov_dir, target_dir, out_dir,
                                  r_targets[i], r_sources[i], r_res[i],
                                  land_mark, theta, exp_d, NA_val, iter,
                                  metrics, connectivity_function_name,
                                  cost_function_name, sensitivity_wrt,
                                  sensitivity_method,
                                  sensitivity_landscape_measure,
                                  sensitivity_unitless,
                                  sensitivity_one_out_of,
                                  sensitivity_diagvalue,
                                  sensitivity_target_equal_source,
                                  sensitivity_require_landmark_one)

                if result == "done"
                    success = true
                    results[i] = "ok"

                    if progress
                        lock(print_lock) do
                            progress_counter[] += 1
                            pct = round(Int, 100 * progress_counter[] / n)
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
                @warn "Failed tile $(r_targets[i]) (index $i) after $max_retries attempt(s): $last_error"
            end
            results[i] = last_error
        end
    end

    return results
end
