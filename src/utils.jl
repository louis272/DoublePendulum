include("system.jl")

"""
    Calcule les dérivées temporelles de l'état (vitesses et accélérations).

    Arguments:
        - dp : L'objet pendule (pour m1, m2, l1, l2, g)
        - state : L'état [θ1, θ2, ω1, ω2] pour lequel on veut la dérivée.

    Retourne:
        - Le vecteur dérivé [dθ1/dt, dθ2/dt, dω1/dt, dω2/dt] = [ω1, ω2, α1, α2]
"""
function derivees(dp::DoublePendule, state::Vector{Float64})
    theta1, theta2, omega1, omega2 = state
    m1 = dp.p1.m
    m2 = dp.p2.m
    l1 = dp.p1.l
    l2 = dp.p2.l
    g  = dp.g

    # Calculs préliminaires
    s1, c1 = sin(theta1), cos(theta1)
    s2, c2 = sin(theta2), cos(theta2)
    delta_theta = theta1 - theta2
    s_delta = sin(delta_theta)
    c_delta = cos(delta_theta)

    denom = 2*m1 + m2 - m2 * cos(2*delta_theta) # Dénominateur commun

    # Calcul des accélérations

    # Accélération angulaire 1
    # Numérateur (Force: Gravité + Centrifuge + Coriolis)
    num_1 = (-g * (2*m1 + m2) * s1
                - m2 * g * sin(theta1 - 2*theta2)
                - 2 * s_delta * m2 * (omega2^2 * l2 + omega1^2 * l1 * c_delta))

    alpha1 = num_1 / (l1 * denom)

    # Accélération angulaire 2
    num_2 = (2 * s_delta * (omega1^2 * l1 * (m1 + m2)
                + g * (m1 + m2) * c1
                + omega2^2 * l2 * m2 * c_delta))

    alpha2 = num_2 / (l2 * denom)

    return [omega1, omega2, alpha1, alpha2]
end


"""
    Effectue un pas d'intégration temporelle (RK4) et met à jour l'état du pendule.

    Arguments:
        - dp : L'objet mutable DoublePendule.
        - dt : Pas de temps.
"""
function rk4_step!(dp::DoublePendule, dt::Float64)
    curr_state = dp.state

    # Calcul des coefficients RK4
    k1 = derivees(dp, curr_state)                 # Pente au début de l'intervalle
    k2 = derivees(dp, curr_state .+ 0.5dt .* k1)  # Pente au milieu de l'intervalle (estimation 1)
    k3 = derivees(dp, curr_state .+ 0.5dt .* k2)  # Pente au milieu de l'intervalle (estimation 2)
    k4 = derivees(dp, curr_state .+ dt .* k3)     # Pente à la fin de l'intervalle

    # Moyenne pondérée des pentes
    dp.state = curr_state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
end

