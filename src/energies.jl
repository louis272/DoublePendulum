include("system.jl")

function kinetic_energy(dp::DoublePendulum)
    """
    Calculates the kinetic energy (T) of the system.
    Formula: T = 0.5*(m1+m2)*l1^2*omega1^2 + 0.5*m2*l2^2*omega2^2 + m2*l1*l2*omega1*omega2*cos(theta1-theta2)

    Args:
        dp: A DoublePendulum object.

    Returns:
        The total kinetic energy T.
    """

    l1 = dp.p1.l
    l2 = dp.p2.l
    m1 = dp.p1.m
    m2 = dp.p2.m
    theta1, theta2, omega1, omega2 = dp.state
    g = dp.g

    ke = (0.5 * (m1 + m2) * l1^2 * omega1^2
            + 0.5 * m2 * l2^2 * omega2^2
            + m2 * l1 * l2 * omega1 * omega2 * cos(theta1 - theta2))

    return ke
end

function potential_energy(dp::DoublePendulum)
    """
    Calculates the potential energy (P) of the system.
    Formula: V = m1*g*y1 + m2*g*y2 with y1 = -l1*cos(theta1) and y2 = y1 - l2*cos(theta2)

    Args:
        dp: A DoublePendulum object.

    Returns:
        The total potential energy P.
    """

    l1 = dp.p1.l
    l2 = dp.p2.l
    m1 = dp.p1.m
    m2 = dp.p2.m
    theta1, theta2, omega1, omega2 = dp.state
    g = dp.g

    y1 = -l1 * cos(theta1)
    y2 = y1 - l2 * cos(theta2)

    return m1 * g * y1 + m2 * g * y2
end

function total_energy(dp::DoublePendulum)
    """
    Returns the total mechanical energy E = T + P.

    Args:
        dp: A DoublePendulum object.

    Returns:
        The total mechanical energy E.
    """
    return kinetic_energy(dp) + potential_energy(dp)
end