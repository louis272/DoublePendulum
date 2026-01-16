mutable struct Pendulum
    l::Float64  # Length [m]
    m::Float64  # Mass [kg]
end

mutable struct DoublePendulum
    p1::Pendulum
    p2::Pendulum
    g::Float64             # Gravitational acceleration [m/s^2]
    state::Vector{Float64} # [θ1, θ2, ω1, ω2]
end

function create_double_pendulum(l1, l2, m1, m2, g, state_init)
    """
    Create and return a DoublePendulum object.

    Args:
        l1 (Float64): Length of the first pendulum rod [m].
        l2 (Float64): Length of the second pendulum rod [m].
        m1 (Float64): Mass of the first pendulum [kg].
        m2 (Float64): Mass of the second pendulum [kg].
        g (Float64): Gravitational acceleration [m/s^2].
        state_init (Vector{Float64}): Initial state [θ1, θ2, ω1, ω2].

    Returns:
        DoublePendulum: A new DoublePendulum object.
    """

    p1 = Pendulum(l1, m1)
    p2 = Pendulum(l2, m2)

    return DoublePendulum(p1, p2, g, state_init)
end

function create_real_double_pendulum()
    """
    Create a DoublePendulum object similar to the one in the video.

    Returns:
        DoublePendulum: A DoublePendulum object with predefined parameters.
    """

    l1 = 91.74e-3 # Length of the first rod [m]
    l2 = 69.33e-3 # Length of the second rod [m]
    #l2 = 67.95e-3 # Adjusted length of the second rod [m]

    m1 = 30.00e-3 # Mass of the first pendulum [kg]
    m2 = 2.00e-3  # Mass of the second pendulum [kg]

    g = 9.81 # Gravitational acceleration [m/s^2]

    theta1_0 = -3.111 # Initial angle 1 [rad]
    theta2_0 = -3.031 # Initial angle 2 [rad]
    omega1_0 = 0.302 # Initial angular velocity 1 [rad/s]
    omega2_0 = 3.567 # Initial angular velocity 2 [rad/s]
    initial_state = [theta1_0, theta2_0, omega1_0, omega2_0]

    return create_double_pendulum(l1, l2, m1, m2, g, initial_state)
end