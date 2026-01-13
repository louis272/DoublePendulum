using Plots
using DelimitedFiles

include("utils.jl")
include("energies.jl")

function mode_live()
    dp = create_real_double_pendulum()

    # Display parameters
    dt = 0.0005          # Calculation step [s]
    steps_per_frame = 40 # Perform 40 calculations before drawing a frame
    L = (dp.p1.l + dp.p2.l) * 1.1 # Graphical limits

    dt_render = steps_per_frame * dt  # Physical time per frame
    sim_time = 0.0 # Total simulation time [s]

    println("Starting live simulation (CTRL+C to stop)")
    try
        while true
            start_time = time()

            # Physical calculations
            for _ in 1:steps_per_frame
                rk4_step!(dp, dt)

                sim_time += dt
            end

            # Retrieve coordinates
            x1, y1, x2, y2 = polar_to_cartesian(dp)

            # Display
            p = plot([0, x1, x2], [0, y1, y2],
                xlims=(-L, L), ylims=(-L, L), aspect_ratio=:equal,
                lw=3, marker=:circle, markersize=8, color=:black,
                legend=false, grid=false, axis=false,
                title="Double Pendulum (Live)"
            )

            time_str = "t = $(round(sim_time, digits=2)) s"
            #annotate!(p, -L*0.8, L*0.9, text(time_str, :blue, 12, :left))

            display(p)

            elapsed_time = time() - start_time
            sleep_duration = max(0.0, dt_render - elapsed_time)
            println(sleep_duration)

            sleep(sleep_duration)
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

    dt = 0.0001 # Time step [s]
    t_max = 10.0 # Simulation duration [s]
    n_steps = Int(floor(t_max / dt))

    ### Store results
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

    # Plot angles over time
    p1 = plot(times, [history_theta1, history_theta2],
        label=["Theta 1" "Theta 2"],
        xlabel="Time [s]", ylabel="Angle [rad]",
        title="Temporal Evolution"
    )

    # Plot and save the variation of total energy over time
    e0 = energies[1]
    energy_deviation = (energies .- e0) ./ abs(e0)

    # Columns: Time | Energy [J] | Relative Error
    data_to_save = [times energies energy_deviation]
    header = ["time_s" "energy_joules" "relative_error"]

    open("./res/energy_stability.csv", "w") do io
        writedlm(io, header, ',')
        writedlm(io, data_to_save, ',')
    end
    println("CSV file successfully generated.")

    p2 = plot(times, energy_deviation,
        label="Relative Error",
        xlabel="Time [s]", ylabel="ΔE / E0",
        title="Energy Stability (Precision)",
        color=:red, lw=1.5
    )

    savefig(p1, "./res/angles_evolution.png")
    println("Graph saved: angles_evolution.png")
    savefig(p2, "./res/energy_stability.png")
    println("Graph saved: energy_stability.png")

    return times, results
end

function main()
    println("Choose the operating mode:")
    println("1. Real-time simulation")
    println("2. Analysis and result saving")
    print("Enter 1 or 2: ")
    choix = readline()

    if choix == "1"
        mode_live()
    elseif choix == "2"
        mode_analyse()
    else
        println("Invalid choice. Please restart the program.")
    end
end

main()