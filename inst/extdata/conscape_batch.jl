using Base.Threads
using LinearAlgebra
using Logging
using ConScape
using SparseArrays
using Statistics

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

function conscape_batch(src_dir, mov_dir, target_dir, out_dir,
                        r_targets::Vector{String}, r_sources::Vector{String},
                        r_res::Vector{String}, land_mark, theta, exp_d, NA_val,
                        max_retries::Int, progress::Bool,
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
                        sensitivity_require_landmark_one = true)

    _set_blas_threads(blas_threads)

    # Align filenames across target/source/movement by basename.
    # This prevents accidental mixing of tiles when some files have been deleted.
    bn = sort!(collect(intersect(intersect(Set(r_targets), Set(r_sources)), Set(r_res))))
    r_targets = bn
    r_sources = bn
    r_res     = bn

    n = length(r_targets)
    results = Vector{Union{String, Nothing}}(undef, n)
    diagnostics = Vector{NamedTuple}(undef, n)
    print_lock = ReentrantLock()
    progress_counter = Ref(0)
    work_estimates = [
        filesize(joinpath(target_dir, r_targets[i])) +
        filesize(joinpath(src_dir, r_sources[i])) +
        filesize(joinpath(mov_dir, r_res[i]))
        for i in 1:n
    ]
    job_order = sortperm(work_estimates, rev = true)
    next_job = Atomic{Int}(1)

    @threads for worker_id in 1:Threads.nthreads()
        while true
            pos = atomic_add!(next_job, 1)
            pos > n && break
            i = job_order[pos]

        # Use the filename stem (e.g., r_191) for a stable, interpretable iteration tag.
        iter = "-" * splitext(r_targets[i])[1]
        attempt = 1
        success = false
        last_error = ""
        started = time()

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
                                  sensitivity_require_landmark_one,
                                  false)

                if result in ("done", "skipped_empty_targets")
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

        elapsed = round(time() - started; digits = 3)
        if !success
            lock(print_lock) do
                @warn "Failed tile $(r_targets[i]) (index $i) after $max_retries attempt(s): $last_error"
            end
            results[i] = last_error
        end

        diagnostics[i] = (
            index = i,
            tile = r_targets[i],
            status = success ? "ok" : "failed",
            attempts = attempt - 1,
            seconds = elapsed,
            work_estimate = work_estimates[i],
            error = success ? "" : last_error
        )
        end
    end

    _write_batch_diagnostics(out_dir, diagnostics)
    return results
end
