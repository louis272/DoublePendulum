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
        history_theta1[i] = dp.state[1] #% (2*pi)
        history_theta2[i] = dp.state[2] #% (2*pi)

        energies[i] = total_energy(dp)
    end

    println("Analysis simulation completed")

    # Plot angles over time
    p1 = plot(times, [wrap_angles(history_theta1), wrap_angles(history_theta2)],
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

function mode_comparaison()
    # Load data from CSV files
    data_path = "./res/video_data.csv"
    if !isfile(data_path)
        println("Error: data file not found.")
        return
    end

    # Read the CSV file
    data = readdlm(data_path, ';', skipstart=1)

    t_exp = data[:, 1]       # Time from experiment
    theta1_exp = data[:, 2]  # Theta 1 from experiment
    theta2_exp = data[:, 3]  # Theta 2 from experiment

    # Unwrap experimental angles
    theta1_exp_unwrapped = zeros(length(theta1_exp))
    theta2_exp_unwrapped = zeros(length(theta2_exp))
    theta1_exp_unwrapped[1] = theta1_exp[1]
    theta2_exp_unwrapped[1] = theta2_exp[1]

    for i in 2:length(t_exp)
        theta1_exp_unwrapped[i] = unwrap_angle(theta1_exp[i], theta1_exp_unwrapped[i-1])
        theta2_exp_unwrapped[i] = unwrap_angle(theta2_exp[i], theta2_exp_unwrapped[i-1])
    end

    # Initialize the simulation
    theta1_0 = theta1_exp_unwrapped[1] # Initial angle theta 1 [rad]
    theta2_0 = theta2_exp_unwrapped[1] # Initial angle theta 2 [rad]
    omega1_0 = (theta1_exp_unwrapped[2] - theta1_exp_unwrapped[1]) / (t_exp[2] - t_exp[1])  # Initial angular velocity omega 1 [rad/s]
    omega2_0 = (theta2_exp_unwrapped[2] - theta2_exp_unwrapped[1]) / (t_exp[2] - t_exp[1])  # Initial angular velocity omega 2 [rad/s]

    # Calculate initial velocities by averaging over N points to avoid noise
    N = min(10, length(t_exp) - 1)  # Use 10 points or fewer if not enough data

    # Angular velocity 1 (average of derivatives)
    omega1_sum = 0.0
    for i in 1:N
        omega1_sum += (theta1_exp_unwrapped[i+1] - theta1_exp_unwrapped[i]) / (t_exp[i+1] - t_exp[i])
    end
    omega1_0 = omega1_sum / N

    # Angular velocity 2 (average of derivatives)
    omega2_sum = 0.0
    for i in 1:N
        omega2_sum += (theta2_exp_unwrapped[i+1] - theta2_exp_unwrapped[i]) / (t_exp[i+1] - t_exp[i])
    end
    omega2_0 = omega2_sum / N

    initial_state = [theta1_0, theta2_0, omega1_0, omega2_0]
    dp = create_real_double_pendulum()
    dp.state = initial_state

    # Simulation
    dt = 0.0001 # Time step [s]
    t_max = maximum(t_exp) # Simulation duration [s]
    n_steps = Int(floor(t_max / dt))

    sim_time = zeros(n_steps)
    sim_theta1 = zeros(n_steps)
    sim_theta2 = zeros(n_steps)

    for i in 1:n_steps
        rk4_step!(dp, dt)
        sim_time[i] = (i) * dt
        sim_theta1[i] = dp.state[1] #% (2*pi)
        sim_theta2[i] = dp.state[2] #% (2*pi)
    end

    # Plot comparison
    # Angle 1
    p1 = plot(t_exp, theta1_exp, seriestype=:scatter, label="Video", markersize=2, color=:blue, alpha=0.5)
    plot!(p1, sim_time, wrap_angles(sim_theta1), label="Simulation", lw=2, color=:red)
    title!(p1, "Theta 1 : Reality vs Simulation")
    ylabel!(p1, "Angle [rad]")

    # Angle 2
    p2 = plot(t_exp, theta2_exp, seriestype=:scatter, label="Video", markersize=2, color=:blue, alpha=0.5)
    plot!(p2, sim_time, wrap_angles(sim_theta2), label="Simulation", lw=2, color=:red)
    title!(p2, "Theta 2 : Reality vs Simulation")
    xlabel!(p2, "Time [s]")
    ylabel!(p2, "Angle [rad]")

    final_plot = plot(p1, p2, layout=(2,1), size=(800,600))
    savefig(final_plot, "./res/comparison_reality_simulation.png")
    println("Graph saved: comparison_reality_simulation.png")
end

function mode_comparaison_2()
    # Load data from CSV files
    data_path = "./res/video_data.csv"
    if !isfile(data_path)
        println("Error: data file not found.")
        return
    end

    # Read the CSV file
    data = readdlm(data_path, ';', skipstart=1)

    t_exp = data[:, 1]       # Time from experiment
    theta1_exp = data[:, 2]  # Theta 1 from experiment (wrapped)
    theta2_exp = data[:, 3]  # Theta 2 from experiment (wrapped)

    # Unwrap experimental angles
    theta1_exp_unwrapped = zeros(length(theta1_exp))
    theta2_exp_unwrapped = zeros(length(theta2_exp))
    theta1_exp_unwrapped[1] = theta1_exp[1]
    theta2_exp_unwrapped[1] = theta2_exp[1]

    for i in 2:length(t_exp)
        theta1_exp_unwrapped[i] = unwrap_angle(theta1_exp[i], theta1_exp_unwrapped[i-1])
        theta2_exp_unwrapped[i] = unwrap_angle(theta2_exp[i], theta2_exp_unwrapped[i-1])
    end

    # Calculate initial conditions from unwrapped data
    theta1_0 = theta1_exp_unwrapped[1]
    theta2_0 = theta2_exp_unwrapped[1]

    # Calculate initial velocities by averaging over N points to reduce noise
    N = min(10, length(t_exp) - 1)

    omega1_sum = 0.0
    for i in 1:N
        omega1_sum += (theta1_exp_unwrapped[i+1] - theta1_exp_unwrapped[i]) / (t_exp[i+1] - t_exp[i])
    end
    omega1_0 = omega1_sum / N

    omega2_sum = 0.0
    for i in 1:N
        omega2_sum += (theta2_exp_unwrapped[i+1] - theta2_exp_unwrapped[i]) / (t_exp[i+1] - t_exp[i])
    end
    omega2_0 = omega2_sum / N

    println("\n=== Initial Conditions ===")
    println("θ1(0) = $(round(theta1_0, digits=3)) rad = $(round(rad2deg(theta1_0), digits=1))°")
    println("θ2(0) = $(round(theta2_0, digits=3)) rad = $(round(rad2deg(theta2_0), digits=1))°")
    println("ω1(0) = $(round(omega1_0, digits=3)) rad/s")
    println("ω2(0) = $(round(omega2_0, digits=3)) rad/s")

    # Initialize simulation
    initial_state = [theta1_0, theta2_0, omega1_0, omega2_0]
    dp = create_real_double_pendulum()
    dp.state = initial_state

    # Run simulation
    dt = 0.0001
    t_max = maximum(t_exp)
    n_steps = Int(floor(t_max / dt))

    sim_time = zeros(n_steps)
    sim_theta1 = zeros(n_steps)
    sim_theta2 = zeros(n_steps)

    println("Running simulation for $(round(t_max, digits=2))s...")

    for i in 1:n_steps
        rk4_step!(dp, dt)
        sim_time[i] = i * dt
        sim_theta1[i] = dp.state[1]  # Keep unwrapped
        sim_theta2[i] = dp.state[2]  # Keep unwrapped
    end

    println("Simulation completed")

    # Plot comparison (all unwrapped for continuous curves)

    # Theta 1
    p1 = plot(t_exp, theta1_exp_unwrapped,
        seriestype=:scatter, label="Experimental Data",
        markersize=2, color=:blue, alpha=0.5)
    plot!(p1, sim_time, sim_theta1,
        label="Simulation", lw=2, color=:red)
    title!(p1, "Theta 1: Reality vs Simulation")
    ylabel!(p1, "Angle [rad]")

    # Theta 2
    p2 = plot(t_exp, theta2_exp_unwrapped,
        seriestype=:scatter, label="Experimental Data",
        markersize=2, color=:blue, alpha=0.5)
    plot!(p2, sim_time, sim_theta2,
        label="Simulation", lw=2, color=:red)
    title!(p2, "Theta 2: Reality vs Simulation")
    xlabel!(p2, "Time [s]")
    ylabel!(p2, "Angle [rad]")

    final_plot = plot(p1, p2, layout=(2,1), size=(800,600))
    savefig(final_plot, "./res/comparison_reality_simulation.png")
    println("Graph saved: comparison_reality_simulation.png")

    # Calculate and plot divergence
    calculate_divergence(t_exp, theta1_exp_unwrapped, theta2_exp_unwrapped,
                        sim_time, sim_theta1, sim_theta2)
end

function calculate_divergence(t_exp, theta1_exp, theta2_exp, sim_time, sim_theta1, sim_theta2)
    """Calculate and plot the divergence between simulation and experiment."""

    distances = zeros(length(sim_time))

    for i in 1:length(sim_time)
        # Find closest experimental time point
        idx_exp = argmin(abs.(t_exp .- sim_time[i]))

        # Euclidean distance in angle space
        distances[i] = sqrt((sim_theta1[i] - theta1_exp[idx_exp])^2 +
                           (sim_theta2[i] - theta2_exp[idx_exp])^2)
    end

    p_div = plot(sim_time, distances,
        label="Divergence",
        xlabel="Time [s]",
        ylabel="Distance [rad]",
        title="Chaotic Divergence",
        lw=2, color=:green,
        legend=:topleft)

    savefig(p_div, "./res/divergence.png")
    println("Graph saved: divergence.png")

    # Estimate Lyapunov exponent (exponential growth phase)
    mask = (sim_time .>= 0.2) .& (sim_time .<= 0.8) .& (distances .> 1e-3)
    if sum(mask) > 10
        log_dist = log.(distances[mask])
        t_subset = sim_time[mask]
        n = length(t_subset)

        # Linear regression: log(distance) = λ*t + const
        lyap = (n * sum(t_subset .* log_dist) - sum(t_subset) * sum(log_dist)) /
               (n * sum(t_subset.^2) - sum(t_subset)^2)

        println("Estimated Lyapunov exponent: λ ≈ $(round(lyap, digits=2)) s⁻¹")
    end
end

function main()
    println("Choose the operating mode:")
    println("1. Real-time simulation")
    println("2. Analysis and result saving")
    println("3. Comparison with experimental data")
    print("Enter 1, 2 or 3 : ")
    choice = readline()

    if choice == "1"
        mode_live()
    elseif choice == "2"
        mode_analyse()
    elseif choice == "3"
        mode_comparaison_2()
    else
        println("Invalid choice. Please restart the program.")
    end
end

main()