using Test
using DecoratedRandomPlanarMaps

@testset "Tutte layouts" begin
    edges = Int32[
        0 1;
        1 2;
        2 3;
        3 0;
        0 4;
        1 4;
        2 4;
        3 4;
    ]
    boundary = Int32[0, 1, 2, 3]
    boundary_pos = Float64[
        -1.0 -1.0;
         1.0 -1.0;
         1.0  1.0;
        -1.0  1.0;
    ]

    pos = compute_tutte_layout(
        5,
        edges,
        boundary;
        boundary_positions=boundary_pos,
        solver="direct",
    )

    @test size(pos) == (5, 2)
    @test pos[Int.(boundary) .+ 1, :] ≈ boundary_pos
    @test pos[5, :] ≈ zeros(2) atol=1.0e-10
end

@testset "Harmonic boundary inference" begin
    edges = Int32[
        0 1;
        1 2;
        2 0;
        0 3;
        1 3;
        2 3;
    ]
    boundary = Int32[0, 1, 2]

    pos, meta = compute_tutte_layout(
        4,
        edges,
        boundary;
        seed=5,
        solver="direct",
        return_metadata=true,
    )

    @test size(pos) == (4, 2)
    @test meta["boundary_positions_mode"] == "harmonic_measure"
    @test haskey(meta, "boundary_info")
end


@testset "AMG solver name is rejected explicitly" begin
    edges = Int32[
        0 1;
        1 2;
        2 0;
    ]
    boundary = Int32[0, 1, 2]
    @test_throws ArgumentError compute_tutte_layout(3, edges, boundary; solver="amg")
end
