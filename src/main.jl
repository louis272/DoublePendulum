using Plots
using DelimitedFiles

include("utils.jl")
include("energies.jl")

function mode_live()
    dp = create_real_double_pendulum()

    # Paramètres d'affichage
    dt = 0.0005          # Pas de calcul
    steps_per_frame = 40 # On fait 40 calculs avant de dessiner une image
    L = (dp.p1.l + dp.p2.l) * 1.1 # Limites graphiques

    println("Starting live simulation (CTRL+C to stop)")
    try
        while true
            # Calculs physiques
            for _ in 1:steps_per_frame
                rk4_step!(dp, dt)
            end

            # Récupération des coordonnées
            x1, y1, x2, y2 = polar_to_cartesian(dp)

            # Affichage
            p = plot([0, x1, x2], [0, y1, y2],
                xlims=(-L, L), ylims=(-L, L), aspect_ratio=:equal,
                lw=3, marker=:circle, markersize=8, color=:black,
                legend=false, grid=false, axis=false,
                title="Double Pendule (Live)"
            )
            display(p)
            sleep(0.005)
        end
    catch e
        if isa(e, InterruptException)
            println("Simulation stopped")
        else
            rethrow(e)
        end
    end
end

function mode_analyse()
    dp = create_real_double_pendulum()

    dt = 0.0001 # Pas de temps (s)
    t_max = 10.0 # Durée simulation (s)
    n_steps = Int(floor(t_max / dt))

    ### Stockage des résultats
    times = zeros(n_steps)
    results = zeros(4, n_steps)
    history_theta1 = zeros(n_steps)
    history_theta2 = zeros(n_steps)
    energies = zeros(n_steps)

    println("Starting analysis simulation for $t_max seconds")

    for i in 1:n_steps
        rk4_step!(dp, dt)

        times[i] = (i) * dt
        results[:, i] = dp.state
        history_theta1[i] = dp.state[1] % (2*pi)
        history_theta2[i] = dp.state[2] % (2*pi)

        energies[i] = total_energy(dp)
    end

    println("Analysis simulation completed")

    # Graphique des angles en fonction du temps
    p1 = plot(times, [history_theta1, history_theta2],
        label=["Theta 1" "Theta 2"],
        xlabel="Temps (s)", ylabel="Angle (rad)",
        title="Evolution temporelle"
    )

    # Graphique et csv de la variation de l'énergie totale en fonction du temps
    e0 = energies[1]
    energy_deviation = (energies .- e0) ./ abs(e0)

    # Colonnes : Temps | Energie (J) | Erreur Relative
    data_to_save = [times energies energy_deviation]
    header = ["temps_s" "energie_joules" "erreur_relative"]

    open("./res/energy_stability.csv", "w") do io
        writedlm(io, header, ',')
        writedlm(io, data_to_save, ',')
    end
    println("Fichier CSV généré avec succès.")

    p2 = plot(times, energy_deviation,
        label="Erreur relative",
        xlabel="Temps (s)", ylabel="ΔE / E0",
        title="Stabilité de l'énergie (Précision)",
        color=:red, lw=1.5
    )

    savefig(p1, "./res/angles_evolution.png")
    println("Graphique sauvegardé : angles_evolution.png")
    savefig(p2, "./res/energy_stability.png")
    println("Graphique sauvegardé : energy_stability.png")

    return times, results
end

function main()
    println("Choisissez le mode de fonctionnement :")
    println("1. Simulation en temps réel")
    println("2. Analyse et sauvegarde des résultats")
    print("Entrez 1 ou 2 : ")
    choix = readline()

    if choix == "1"
        mode_live()
    elseif choix == "2"
        mode_analyse()
    else
        println("Choix invalide. Veuillez relancer le programme.")
    end
end

main()