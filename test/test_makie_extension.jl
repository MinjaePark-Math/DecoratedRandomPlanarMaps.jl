using Test
using DecoratedRandomPlanarMaps

const MAKIE_PIPELINE_VIEWER_CALLS = Ref(0)
const MAKIE_PIPELINE_LAST_RADII = Ref{Any}(nothing)
const MAKIE_PIPELINE_3D_VIEWER_CALLS = Ref(0)

@testset "Makie extension loads on upgraded environment" begin
    loaded = false
    load_error = nothing

    try
        @eval using GLMakie
        @eval using GeometryBasics
        loaded = true
    catch err
        load_error = sprint(showerror, err)
    end

    if !loaded
        @info "Skipping GLMakie optional backend check" julia_version=string(VERSION) error=load_error
        @test load_error !== nothing
    else
        ext = Base.get_extension(DecoratedRandomPlanarMaps, :DecoratedRandomPlanarMapsMakieExt)
        @test ext !== nothing
        @test hasmethod(
            DecoratedRandomPlanarMaps.render_makie_2d,
            Tuple{Any, AbstractMatrix},
        )

        @eval DecoratedRandomPlanarMaps begin
            function render_makie_2d(map_data::SchnyderMap, pos::AbstractMatrix; circle_radii=nothing, kwargs...)
                Main.MAKIE_PIPELINE_VIEWER_CALLS[] += 1
                Main.MAKIE_PIPELINE_LAST_RADII[] = circle_radii
                return nothing
            end

            function render_makie_3d(map_data::SchnyderMap, pos::AbstractMatrix; kwargs...)
                Main.MAKIE_PIPELINE_3D_VIEWER_CALLS[] += 1
                return nothing
            end
        end

        @testset "2D show opens the Makie viewer" begin
            Main.MAKIE_PIPELINE_VIEWER_CALLS[] = 0
            Main.MAKIE_PIPELINE_LAST_RADII[] = nothing

            mktempdir() do tmp
                cd(tmp) do
                    timings = DecoratedRandomPlanarMaps.run_pipeline(Dict(
                        "model" => Dict("type" => "schnyder", "faces" => 8, "seed" => 7),
                        "layout" => Dict("dimension" => 2, "engine" => "circle_packing"),
                        "output" => Dict("show" => true),
                    ))

                    @test Main.MAKIE_PIPELINE_VIEWER_CALLS[] == 1
                    @test Main.MAKIE_PIPELINE_LAST_RADII[] !== nothing
                    @test all(step.name != "preview_svg" for step in timings.steps)
                    @test !isfile(joinpath(tmp, "preview.svg"))
                end
            end
        end

        @testset "3D show opens the Makie viewer" begin
            Main.MAKIE_PIPELINE_3D_VIEWER_CALLS[] = 0

            mktempdir() do tmp
                cd(tmp) do
                    timings = DecoratedRandomPlanarMaps.run_pipeline(Dict(
                        "model" => Dict("type" => "schnyder", "faces" => 8, "seed" => 7),
                        "layout" => Dict("dimension" => 3, "engine" => "sfdp"),
                        "output" => Dict("show" => true),
                    ))

                    @test Main.MAKIE_PIPELINE_3D_VIEWER_CALLS[] == 1
                    @test all(step.name != "export_web" for step in timings.steps)
                    @test !isdir(joinpath(tmp, "decoratedrandomplanarmaps_preview_web"))
                end
            end
        end
    end
end
