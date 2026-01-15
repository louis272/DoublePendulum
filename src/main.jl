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

function mode_analysis()
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

function mode_comparison()
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

    # Calculate initial angular velocities by averaging over N points to avoid noise
    N = min(10, length(t_exp) - 1)  # Use 10 points or fewer if not enough data

    # Angular velocity 1 (average of derivatives) [rad/s]
    omega1_sum = 0.0
    for i in 1:N
        omega1_sum += (theta1_exp_unwrapped[i+1] - theta1_exp_unwrapped[i]) / (t_exp[i+1] - t_exp[i])
    end
    omega1_0 = omega1_sum / N

    # Angular velocity 2 (average of derivatives) [rad/s]
    omega2_sum = 0.0
    for i in 1:N
        omega2_sum += (theta2_exp_unwrapped[i+1] - theta2_exp_unwrapped[i]) / (t_exp[i+1] - t_exp[i])
    end
    omega2_0 = omega2_sum / N

    println("================================")
    println("Initial Conditions")
    println("θ1(0) = $(round(theta1_0, digits=3)) rad = $(round(rad2deg(theta1_0), digits=1))°")
    println("θ2(0) = $(round(theta2_0, digits=3)) rad = $(round(rad2deg(theta2_0), digits=1))°")
    println("ω1(0) = $(round(omega1_0, digits=3)) rad/s")
    println("ω2(0) = $(round(omega2_0, digits=3)) rad/s")
    println("================================")

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
        sim_theta1[i] = dp.state[1]
        sim_theta2[i] = dp.state[2]
    end

    # Plot comparison (wrapped angles)
    # Angle 1
    p1 = plot(t_exp, theta1_exp, label="Video", lw=2, color=:blue)
    plot!(p1, sim_time, wrap_angles(sim_theta1), label="Simulation", lw=2, color=:red)
    title!(p1, "Theta 1 : Reality vs Simulation")
    ylabel!(p1, "Angle [rad]")

    # Angle 2
    p2 = plot(t_exp, theta2_exp, label="Video", lw=2, color=:blue)
    plot!(p2, sim_time, wrap_angles(sim_theta2), label="Simulation", lw=2, color=:red)
    title!(p2, "Theta 2 : Reality vs Simulation")
    xlabel!(p2, "Time [s]")
    ylabel!(p2, "Angle [rad]")

    final_plot = plot(p1, p2, layout=(2,1), size=(800,600))
    savefig(final_plot, "./res/comparison_reality_simulation_wrapped.png")

    # Plot comparison (angles unwrapped for continuous curves)
    # Theta 1
    p1 = plot(t_exp, theta1_exp_unwrapped, seriestype=:scatter, label="Experimental Data", markersize=2, color=:blue, alpha=0.5)
    plot!(p1, sim_time, sim_theta1, label="Simulation", lw=2, color=:red)
    title!(p1, "Theta 1: Reality vs Simulation")
    ylabel!(p1, "Angle [rad]")

    # Theta 2
    p2 = plot(t_exp, theta2_exp_unwrapped, seriestype=:scatter, label="Experimental Data", markersize=2, color=:blue, alpha=0.5)
    plot!(p2, sim_time, sim_theta2, label="Simulation", lw=2, color=:red)
    title!(p2, "Theta 2: Reality vs Simulation")
    xlabel!(p2, "Time [s]")
    ylabel!(p2, "Angle [rad]")

    final_plot = plot(p1, p2, layout=(2,1), size=(800,600))
    savefig(final_plot, "./res/comparison_reality_simulation_unwrapped.png")
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
        mode_analysis()
    elseif choice == "3"
        mode_comparison()
    else
        println("Invalid choice. Please restart the program.")
    end
end

main()