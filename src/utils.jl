include("system.jl")

function derivees(dp::DoublePendulum, state::Vector{Float64})
    """
    Calculate the time derivatives of the state (velocities and accelerations).

    Args:
        dp : The DoublePendulum object (for m1, m2, l1, l2, g).
        state : The state [θ1, θ2, ω1, ω2] for which the derivative is needed.

    Returns:
        The derivative vector [dθ1/dt, dθ2/dt, dω1/dt, dω2/dt] = [ω1, ω2, α1, α2].
    """

    theta1, theta2, omega1, omega2 = state
    m1 = dp.p1.m
    m2 = dp.p2.m
    l1 = dp.p1.l
    l2 = dp.p2.l
    g  = dp.g

    # Preliminary calculations
    s1, c1 = sin(theta1), cos(theta1)
    s2, c2 = sin(theta2), cos(theta2)
    delta_theta = theta1 - theta2
    s_delta = sin(delta_theta)
    c_delta = cos(delta_theta)

    denom = 2*m1 + m2 - m2 * cos(2*delta_theta) # Common denominator

    # Calculate accelerations

    # Angular acceleration 1
    # Numerator (Force: Gravity + Centrifugal + Coriolis)
    num_1 = (-g * (2*m1 + m2) * s1
                - m2 * g * sin(theta1 - 2*theta2)
                - 2 * s_delta * m2 * (omega2^2 * l2 + omega1^2 * l1 * c_delta))

    alpha1 = num_1 / (l1 * denom)

    # Angular acceleration 2
    num_2 = (2 * s_delta * (omega1^2 * l1 * (m1 + m2)
                + g * (m1 + m2) * c1
                + omega2^2 * l2 * m2 * c_delta))

    alpha2 = num_2 / (l2 * denom)

    return [omega1, omega2, alpha1, alpha2]
end

function rk4_step!(dp::DoublePendulum, dt::Float64)
    """
    Perform a time integration step (RK4) and update the state of the DoublePendulum.

    Args:
        dp : The mutable DoublePendulum object.
        dt : Time step.
    """

    curr_state = dp.state

    # Calculate RK4 coefficients
    k1 = derivees(dp, curr_state)                 # Slope at the beginning of the interval
    k2 = derivees(dp, curr_state .+ 0.5dt .* k1)  # Slope in the middle of the interval (estimate 1)
    k3 = derivees(dp, curr_state .+ 0.5dt .* k2)  # Slope in the middle of the interval (estimate 2)
    k4 = derivees(dp, curr_state .+ dt .* k3)     # Slope at the end of the interval

    # Weighted average of slopes
    dp.state = curr_state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
end

function polar_to_cartesian(dp::DoublePendulum)
    """
    Convert polar coordinates to cartesian coordinates for visualization.

    Args:
        dp : A DoublePendulum object.

    Returns:
        The coordinates (x1, y1, x2, y2) of the DoublePendulum.
    """

    l1 = dp.p1.l
    l2 = dp.p2.l
    theta1 = dp.state[1]
    theta2 = dp.state[2]

    x1 = l1 * sin(theta1)
    y1 = -l1 * cos(theta1)

    x2 = x1 + l2 * sin(theta2)
    y2 = y1 - l2 * cos(theta2)

    return (x1, y1, x2, y2)
end

function wrap_angle(angle::Float64)
    """
    Wrap an angle to the range [-π, π].

    Args:
        angle : The angle in radians.

    Returns:
        The wrapped angle in radians.
    """

    return mod(angle + π, 2π) - π
end

function wrap_angles(angles::Vector{Float64})
    """
    Wrap a vector of angles to the range [-π, π].

    Args:
        angles : A vector of angles in radians.

    Returns:
        A vector of wrapped angles in radians.
    """

    return [wrap_angle(angle) for angle in angles]
end

function unwrap_angle(angle::Float64, angle_prev::Float64)
    """
    Unwrap an angle, given the previous unwrapped angle.

    Args:
        angle : The angle in radians.

    Returns:
        The unwrapped angle in radians.
    """

    # Calculer la différence entre l'angle actuel et le précédent
    delta = angle - angle_prev

    # Ramener cette différence dans [-π, π]
    delta = mod(delta + pi, 2*pi) - pi

    # L'angle unwrappé est l'angle précédent + la différence corrigée
    return angle_prev + delta
end