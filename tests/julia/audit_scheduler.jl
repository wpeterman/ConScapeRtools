using Test, Distributed, LinearAlgebra, Logging
pool = addprocs(2; exeflags = "--threads=1")
try
    # A yielding stand-in exposes overlapping calls that real CPU-bound solves
    # can conceal. Neither input matrices nor full manifests should be queued.
    @everywhere begin
        using LinearAlgebra
        active_jobs = 0
        peak_jobs = 0
        function _as_string_vector(x)
            x isa AbstractString ? [String(x)] : String.(x)
        end
        function conscape(src, mov, target, out, t::AbstractString,
                          s::AbstractString, m::AbstractString, args...)
            global active_jobs += 1
            global peak_jobs = max(peak_jobs, active_jobs)
            try
                sleep(0.03)
                t == "tile-3.asc" && error("intentional tile failure")
                "done"
            finally
                global active_jobs -= 1
            end
        end
    end
    original_logger = global_logger()
    helper_file = joinpath(abspath(ARGS[1]), "inst", "extdata", "conscape_batch_distributed.jl")
    @everywhere [1; pool] include($helper_file)
    @test global_logger() === original_logger
    mktempdir() do dir
        files = ["tile-$i.asc" for i in 1:12]
        for f in files
            write(joinpath(dir, f), "test")
        end
        before = [remotecall_fetch(BLAS.get_num_threads, p) for p in pool]
        result = conscape_batch_distributed(dir, dir, dir, dir, files, files,
            files, 1, 0.1, 3, 1e-50, 1, false, 1,
            ["connected_habitat"], "expected_cost", "minuslog", String[],
            "analytical", "sum", true, 1, nothing, true, true, pool)
        @test count(==("ok"), result) == 11
        @test all(remotecall_fetch(Core.eval, p, Main, :peak_jobs) == 1 for p in pool)
        @test [remotecall_fetch(BLAS.get_num_threads, p) for p in pool] == before
        @test countlines(joinpath(dir, "conscape_batch_diagnostics.csv")) >= 13
        @test_throws ArgumentError conscape_batch_distributed(dir, dir, dir, dir,
            files, files[1:3], files, 1, 0.1, 3, 1e-50)
    end
finally
    rmprocs(pool)
end
