mutable struct Pendule
    l::Float64  # Longueur [m]
    m::Float64  # Masse [kg]
end

mutable struct DoublePendule
    p1::Pendule
    p2::Pendule
    g::Float64             # Gravité [m/s^2]
    state::Vector{Float64} # [θ1, θ2, ω1, ω2]
end

"""
    Crée et retourne un objet DoublePendule.
"""
function create_double_pendulum(l1, l2, m1, m2, g, state_init)
    p1 = Pendule(l1, m1)
    p2 = Pendule(l2, m2)

    return DoublePendule(p1, p2, g, state_init)
end

"""
    Crée un pendule double similaire à celui de la vidéo.
"""
function create_real_double_pendulum()
    l1 = 91.74e-2 # Longueur tige 1 (m)
    l2 = 69.33e-2 # Longueur tige 2 (m)

    m1 = 30.00e-2 # Masse aléatoire pendule 1 (kg)
    m2 = 2.00e-2  # Masse aléatoire pendule 2 (kg)

    g = 9.81 # Accélération gravitationnelle (m/s²)

    theta1_0 = π/2
    theta2_0 = π/2
    omega1_0 = 0.0
    omega2_0 = 0.0
    initial_state = [theta1_0, theta2_0, omega1_0, omega2_0]

    return create_double_pendulum(l1, l2, m1, m2, g, initial_state)
end