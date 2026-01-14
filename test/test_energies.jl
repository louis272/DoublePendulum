using Test
include("../src/energies.jl")

@testset "Energy Tests" begin

    # Double pendulum parameters for tests
    l1, l2 = 1.0, 1.0
    m1, m2 = 1.0, 1.0
    g = 10.0

    # Helper to create a pendulum
    function make_dp(theta1, theta2, omega1, omega2)
        state = [theta1, theta2, omega1, omega2]
        return create_double_pendulum(l1, l2, m1, m2, g, state)
    end

    @testset "Kinetic Energy (T)" begin
        # Case 1: Zero velocity -> Zero kinetic energy
        dp_rest = make_dp(0.0, 0.0, 0.0, 0.0)
        @test kinetic_energy(dp_rest) == 0.0

        # Case 2: Unit velocities, pendulum aligned (theta1 = theta2 = 0)
        dp_moving = make_dp(0.0, 0.0, 1.0, 1.0)
        @test isapprox(kinetic_energy(dp_moving), 2.5)

        # Case 3: Perpendicular pendulums (theta1=0, theta2=pi/2) -> cos(pi/2) = 0
        dp_ortho = make_dp(0.0, π/2, 1.0, 1.0)
        @test isapprox(kinetic_energy(dp_ortho), 1.5)
    end

    @testset "Potential Energy (V)" begin
        # Case 1: Horizontal (theta = pi/2) -> y1=0, y2=0 -> V=0
        dp_flat = make_dp(π/2, π/2, 0.0, 0.0)
        @test isapprox(potential_energy(dp_flat), 0.0, atol=1e-12)

        # Case 2: Vertical down (theta = 0) -> y1=-1, y2=-2
        dp_down = make_dp(0.0, 0.0, 0.0, 0.0)
        @test isapprox(potential_energy(dp_down), -30.0)

        # Case 3: Vertical up (theta = pi) -> y1=1, y2=2
        dp_up = make_dp(π, π, 0.0, 0.0)
        @test isapprox(potential_energy(dp_up), 30.0)
    end

    @testset "Total Energy (E)" begin
        # E = T + V
        dp_combo = make_dp(0.0, 0.0, 1.0, 1.0)

        @test isapprox(kinetic_energy(dp_combo), 2.5)
        @test isapprox(potential_energy(dp_combo), -30.0)
        @test isapprox(total_energy(dp_combo), -27.5)
    end
end