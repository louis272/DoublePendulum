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