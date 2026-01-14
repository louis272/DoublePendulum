using Test
include("../src/system.jl")

@testset "System Tests" begin

    @testset "Pendulum Struct" begin
        pendulum = Pendulum(1.0, 2.0)
        @test pendulum.l == 1.0
        @test pendulum.m == 2.0
    end

    @testset "DoublePendulum Struct" begin
        pendulum1 = Pendulum(1.0, 2.0)
        pendulum2 = Pendulum(1.5, 3.0)
        state = [π/4, π/6, 0.5, 0.3]
        dp = DoublePendulum(pendulum1, pendulum2, 9.81, state)

        @test dp.p1 == pendulum1
        @test dp.p2 == pendulum2
        @test dp.g == 9.81
        @test dp.state == state
    end

    @testset "create_double_pendulum" begin
        l1, l2, m1, m2, g = 1.0, 1.5, 2.0, 3.0, 9.81
        state_init = [π/4, π/6, 0.5, 0.3]
        dp = create_double_pendulum(l1, l2, m1, m2, g, state_init)

        @test dp.p1.l == l1
        @test dp.p1.m == m1
        @test dp.p2.l == l2
        @test dp.p2.m == m2
        @test dp.g == g
        @test dp.state == state_init
    end
end
