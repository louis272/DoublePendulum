using Test
include("../src/utils.jl")

@testset "Utility Tests" begin

    pi_f = float(π)

    @testset "Angle Management" begin
        # Test wrap_angle
        @test wrap_angle(0.0) == 0.0
        @test isapprox(wrap_angle(pi_f), -pi_f, atol=1e-12)
        @test isapprox(wrap_angle(-pi_f), -pi_f, atol=1e-12)
        @test isapprox(wrap_angle(3pi_f), -pi_f, atol=1e-12)
        @test isapprox(wrap_angle(-3pi_f), -pi_f, atol=1e-12)
        @test isapprox(wrap_angle(2pi_f + 0.1), 0.1, atol=1e-12)

        # Test wrap_angles (Vectors)
        angles = [0.0, 3pi_f, -3pi_f]
        wrapped = wrap_angles(angles)
        @test wrapped[1] == 0.0
        @test isapprox(wrapped[2], -pi_f, atol=1e-12)
        @test isapprox(wrapped[3], -pi_f, atol=1e-12)

        # Test unwrap_angle
        # Case without jump
        @test isapprox(unwrap_angle(0.2, 0.1), 0.2)

        # Case with jump: goes from π - 0.1 to -π + 0.1
        val_prev = pi_f - 0.1
        val_curr = -pi_f + 0.1
        unwrapped = unwrap_angle(val_curr, val_prev)
        @test isapprox(unwrapped, pi_f + 0.1, atol=1e-12)
    end

    @testset "Polar -> Cartesian Conversion" begin
        # Creation of a dummy pendulum
        l1, l2 = 1.0, 1.0
        m1, m2 = 1.0, 1.0
        g = 9.81

        # Case 1: Straight down (0, 0)
        dp_down = create_double_pendulum(l1, l2, m1, m2, g, [0.0, 0.0, 0.0, 0.0])
        x1, y1, x2, y2 = polar_to_cartesian(dp_down)

        @test isapprox(x1, 0.0, atol=1e-12)
        @test isapprox(y1, -1.0, atol=1e-12)
        @test isapprox(x2, 0.0, atol=1e-12)
        @test isapprox(y2, -2.0, atol=1e-12)

        # Case 2: Straight right (π/2, π/2)
        dp_right = create_double_pendulum(l1, l2, m1, m2, g, [π/2, π/2, 0.0, 0.0])
        x1_r, y1_r, x2_r, y2_r = polar_to_cartesian(dp_right)

        @test isapprox(x1_r, 1.0, atol=1e-12)
        @test isapprox(y1_r, 0.0, atol=1e-12)
        @test isapprox(x2_r, 2.0, atol=1e-12)
        @test isapprox(y2_r, 0.0, atol=1e-12)

        # Case 3: Right angle (π/2, π)
        dp_mix = create_double_pendulum(l1, l2, m1, m2, g, [π/2, π, 0.0, 0.0])
        x1_m, y1_m, x2_m, y2_m = polar_to_cartesian(dp_mix)

        @test isapprox(x1_m, 1.0, atol=1e-12)
        @test isapprox(y1_m, 0.0, atol=1e-12)
        @test isapprox(x2_m, 1.0, atol=1e-12)
        @test isapprox(y2_m, 1.0, atol=1e-12)
    end

    @testset "Derivatives Calculation (Physics)" begin
        # Creation of a dummy pendulum
        l1, l2, m1, m2, g = 1.0, 1.0, 1.0, 1.0, 10.0

        # Stationary state: theta=0, omega=0.
        # Angular acceleration must be zero (stable equilibrium)
        dp_stable = create_double_pendulum(l1, l2, m1, m2, g, [0.0, 0.0, 0.0, 0.0])
        res = derivees(dp_stable, dp_stable.state)

        @test length(res) == 4
        @test res[1] == 0.0 # dtheta1/dt = omega1
        @test res[2] == 0.0 # dtheta2/dt = omega2
        @test isapprox(res[3], 0.0, atol=1e-12) # alpha1
        @test isapprox(res[4], 0.0, atol=1e-12) # alpha2

        # Horizontal state (π/2) without velocity
        # Gravity must create immediate acceleration
        dp_falling = create_double_pendulum(l1, l2, m1, m2, g, [π/2, π/2, 0.0, 0.0])
        res_falling = derivees(dp_falling, dp_falling.state)

        @test res_falling[1] == 0.0 # Initial velocity zero
        @test res_falling[3] < 0.0  # Must start falling (negative acceleration for theta)
    end

    @testset "RK4 Integrator" begin
        # Creation of a dummy pendulum
        l1, l2, m1, m2, g = 1.0, 1.0, 1.0, 1.0, 10.0
        dp = create_double_pendulum(l1, l2, m1, m2, g, [π/2, π/2, 0.0, 0.0])

        initial_state = copy(dp.state)
        dt = 0.01
        rk4_step!(dp, dt)

        # After one time step, the state must have changed
        @test dp.state != initial_state

        # Velocity must have increased (in absolute value)
        @test abs(dp.state[3]) > 0.0

        # Test with dt = 0
        state_before_zero = copy(dp.state)
        rk4_step!(dp, 0.0)
        @test dp.state == state_before_zero
    end
end