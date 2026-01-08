include("system.jl")

"""
    Calcule l'énergie cinétique (T) du système.
    Formule : T = 0.5*(m1+m2)*l1^2*omega1^2 + 0.5*m2*l2^2*omega2^2 + m2*l1*l2*omega1*omega2*cos(theta1-theta2)

    Arguments:
        - dp : Un DoublePendule.

    Retourne:
        - L'énergie cinétique totale T.
"""
function kinetic_energy(dp::DoublePendule)
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

"""
    Calcule l'énergie potentielle (P) du système.
    Formule: V = m1*g*y1 + m2*g*y2 avec y1 = -l1*cos(theta1) et y2 = y1 - l2*cos(theta2)

    Arguments:
        - dp : Un DoublePendule.

    Retourne:
        - L'énergie potentielle totale P.
"""
function potential_energy(dp::DoublePendule)
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

"""
    Retourne l'énergie mécanique totale E = T + P.

    Arguments:
        - dp : Un DoublePendule.

    Retourne:
        - L'énergie mécanique totale E.
"""
function total_energy(dp::DoublePendule)
    return kinetic_energy(dp) + potential_energy(dp)
end