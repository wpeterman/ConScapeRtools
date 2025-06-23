using Base.Threads

function conscape_batch(src_dir, mov_dir, target_dir, out_dir,
                       r_targets::Vector{String}, r_sources::Vector{String},
                       r_res::Vector{String}, land_mark, theta, exp_d, NA_val)

    # Initialize thread-safe output and error tracking
    # results = Vector{Any}(undef, length(r_targets))
    errors = Vector{Union{Nothing,Exception}}(nothing, length(r_targets))

    # Use a lock for any thread-unsafe operations (like writing to console/files)
    print_lock = ReentrantLock()

    results = Vector{Union{String,Nothing}}(undef, length(r_targets))

    Threads.@threads for i in 1:length(r_targets)
        iter = "-iter_$i"
        try
            conscape(src_dir, mov_dir, target_dir, out_dir,
                                 r_targets[i], r_sources[i], r_res[i],
                                 land_mark, theta, exp_d, NA_val, iter)

                                             results[i] = "ok"

            # Optional progress reporting (thread-safe)
            lock(print_lock) do
                println("Completed iteration $i/$(length(r_targets))")
            end
        catch e
            results[i] = sprint(showerror, e)
            lock(print_lock) do
                @warn "Failed at iteration $i" exception=(e, catch_backtrace())
            end
        end
    end

    # Check for any failures and return results
    failed_indices = findall(!isnothing, errors)
    if !isempty(failed_indices)
        @warn "Failed iterations: $failed_indices"
    end

    return results
end
