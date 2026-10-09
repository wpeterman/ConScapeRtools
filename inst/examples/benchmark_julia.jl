# julia --startup-file=no benchmark_julia.jl REPO SIDE DENSITY METHOD BLAS REPEATS BLOCK
# METHOD: separate (upstream), shared (one proximity matrix), streamed (blocks).
# Each invocation runs one method, allowing external process-peak measurement.
# Timings include graph construction and all solves, but exclude loading, small
# warm-up, input generation, and ASCII IO. Allocations are not peak RAM.
using ConScape, SparseArrays, LinearAlgebra, Random, Statistics, SHA, Pkg
include(joinpath(abspath(ARGS[1]), "inst", "extdata", "conscape.jl"))

function inputs(n, density)
    rng = MersenneTwister(2718)
    mov = 0.2 .+ 0.6 .* rand(rng, n, n)
    source, target = rand(rng, n, n), rand(rng, n, n)
    target[rand(rng, n, n) .> density] .= 0
    mov, source, sparse(target)
end

function calculate(mov, source, target, method, blocksize)
    g = ConScape.Grid(size(mov)...;
        affinities = ConScape.graph_matrix_from_raster(mov),
        source_qualities = source, target_qualities = target, costs = ConScape.MinusLog())
    transform = x -> exp(-x / 10)
    if method == "streamed"
        return _streamed_default_surfaces(g, 0.1, transform; blocksize)
    end
    h = ConScape.GridRSP(g; θ = 0.1)
    if method == "shared"
        return _paired_metric_surfaces(h, ConScape.expected_cost, transform)
    elseif method == "separate"
        return ConScape.connected_habitat(h; distance_transformation = transform),
               ConScape.betweenness_kweighted(h; distance_transformation = transform)
    end
    error("Unknown method: $method")
end

function main(args)
    n = length(args) >= 2 ? parse(Int, args[2]) : 50
    density = length(args) >= 3 ? parse(Float64, args[3]) : 1.0
    method = length(args) >= 4 ? args[4] : "streamed"
    blas = length(args) >= 5 ? parse(Int, args[5]) : 1
    repeats = length(args) >= 6 ? parse(Int, args[6]) : 3
    blocksize = length(args) >= 7 ? parse(Int, args[7]) : 128
    BLAS.set_num_threads(blas)
    calculate(inputs(8, 1.0)..., method, blocksize)
    mov, source, target = inputs(n, density)
    info = only(p for p in values(Pkg.dependencies()) if p.name == "ConScape")
    println("# Julia=", VERSION, ";ConScape=", pkgversion(ConScape),
            ";tree=", info.tree_hash, ";CPU=", Sys.CPU_NAME,
            ";threads=", Threads.nthreads(), ";BLAS=", BLAS.get_config())
    println("# helper_sha256=", bytes2hex(sha256(read(joinpath(args[1], "inst", "extdata", "conscape.jl")))))
    println("method,side,targets,blas,block,repeat,seconds,allocated_mb,fcon_sum,btwn_sum")
    for rep in 1:repeats
        GC.gc()
        run = @timed calculate(mov, source, target, method, blocksize)
        println(join((method, n, nnz(target), blas, blocksize, rep, run.time,
                      run.bytes / 1e6, sum(run.value[1]), sum(run.value[2])), ','))
    end
end
main(ARGS)
