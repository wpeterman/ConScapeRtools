using Test, ConScape, SparseArrays, LinearAlgebra, Random
include(joinpath(abspath(ARGS[1]), "inst", "extdata", "conscape.jl"))
BLAS.set_num_threads(1)
Random.seed!(2718)

@testset "shared proximity preserves ConScape metric kernels" begin
    for n in (4, 7), sparse_targets in (false, true)
        mov = 0.2 .+ 0.6 .* rand(n, n + 1)
        source = rand(n, n + 1)
        target = rand(n, n + 1)
        sparse_targets && (target[1:2:end] .= 0)
        g = ConScape.Grid(size(mov)...,
            affinities = ConScape.graph_matrix_from_raster(mov),
            source_qualities = source, target_qualities = sparse(target),
            costs = ConScape.MinusLog())
        @test ConScape.coarse_graining(g, 1) == g.target_qualities
        h = ConScape.GridRSP(g; θ = 0.1)
        for conn in (ConScape.expected_cost, ConScape.free_energy_distance,
                     ConScape.least_cost_distance, ConScape.survival_probability,
                     ConScape.power_mean_proximity)
            transform = x -> exp(-x / 3)
            fcon, btwn = _paired_metric_surfaces(h, conn, transform)
            @test fcon ≈ ConScape.connected_habitat(h;
                connectivity_function = conn, distance_transformation = transform)
            @test btwn ≈ ConScape.betweenness_kweighted(h;
                connectivity_function = conn, distance_transformation = transform)
        end
    end
end

@testset "target-block solves preserve full-graph default metrics" begin
    for theta in (0.01, 0.2), costs in (ConScape.MinusLog(), ConScape.Inv()),
        masked in (false, true), dense in (false, true)
        mov = 0.2 .+ 0.6 .* rand(7, 9)
        masked && (mov[1:2, 1:2] .= 0)
        source, target = 3 .* rand(7, 9), 4 .* rand(7, 9)
        target[1:3:end] .= 0
        g = ConScape.Grid(size(mov)...;
            affinities = ConScape.graph_matrix_from_raster(mov),
            source_qualities = source, target_qualities = dense ? target : sparse(target),
            costs = costs)
        h = ConScape.GridRSP(g; θ = theta)
        transform = x -> exp(-x / 3)
        ref = _paired_metric_surfaces(h, ConScape.expected_cost, transform)
        retained = g.id_to_grid_coordinate_list
        for blocksize in (1, 3, 128)
            actual = _streamed_default_surfaces(g, theta, transform; blocksize)
            for i in 1:2
                @test actual[i][retained] ≈ ref[i][retained] rtol=1e-10 atol=1e-10
                @test isnan.(actual[i]) == isnan.(ref[i])
            end
        end
        fcon, btwn = _streamed_default_surfaces(g, theta, transform; between = false)
        @test fcon[retained] ≈ ref[1][retained]
        @test btwn === nothing
        fcon, btwn = _streamed_default_surfaces(g, theta, transform; connected = false)
        @test fcon === nothing
        @test btwn[retained] ≈ ref[2][retained]
    end
end

@testset "coarse graining reports partial-block target loss" begin
    mov = fill(0.5, 5, 7)
    quality = ones(size(mov))
    g = ConScape.Grid(size(mov)...; affinities = ConScape.graph_matrix_from_raster(mov),
        source_qualities = quality, target_qualities = sparse(quality))
    @test_logs (:warn, r"Coarse graining changed") _check_coarse_target_mass(
        g, ConScape.coarse_graining(g, 3), 3)
    @test_throws ArgumentError _check_coarse_target_mass(g, spzeros(size(mov)...), 3)
    @test_logs _check_coarse_target_mass(g, g.target_qualities, 1)
end

@testset "disconnected habitat cannot be silently discarded" begin
    mov = fill(0.5, 5, 7)
    mov[:, 4] .= 0
    quality = ones(size(mov))
    g = ConScape.Grid(size(mov)...; affinities = ConScape.graph_matrix_from_raster(mov),
        source_qualities = quality, target_qualities = sparse(quality))
    @test_throws ArgumentError _check_retained_habitat(g, mov, quality, quality)
end
