using Plots
using DelimitedFiles
using Statistics

include("utils.jl")
include("energies.jl")

function mode_live()
    """
    Run a real-time simulation of the double pendulum, displaying its motion live.
    """

    dp = create_real_double_pendulum()

    # Parameters
    dt = 0.0005                    # Calculation step [s]
    steps_per_frame = 40           # Perform 40 calculations before drawing a frame
    L = (dp.p1.l + dp.p2.l) * 1.1  # Graphical limits

    dt_render = steps_per_frame * dt  # Physical time per frame [s]
    sim_time = 0.0                    # Total simulation time [s]

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
            #annotate!(p, -L*0.8, L*0.9, text(time_str, :blue, 12, :left))  # Decomment to show time on plot

            display(p)

            elapsed_time = time() - start_time
            sleep_duration = max(0.0, dt_render - elapsed_time)

            sleep(sleep_duration)  # Control frame rate
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
    """
    Perform a detailed analysis simulation of the double pendulum over a fixed duration.

    This function calculates and saves the temporal evolution of angles and energy stability to files.

    Returns:
        A tuple containing the simulation times and results.
    """

    dp = create_real_double_pendulum()

    dt = 0.0001                       # Time step [s]
    t_max = 5.0                       # Simulation duration [s]
    n_steps = Int(floor(t_max / dt))  # Number of steps

    ### Store results
    times = zeros(n_steps)           # To store time values [s]
    results = zeros(4, n_steps)      # To store [theta1, theta2, omega1, omega2] over time
    history_theta1 = zeros(n_steps)  # To store theta 1 over time [rad]
    history_theta2 = zeros(n_steps)  # To store theta 2 over time [rad]
    energies = zeros(n_steps)        # To store total energy over time [J]

    # Store trajectory for animation
    traj_x = zeros(n_steps)
    traj_y = zeros(n_steps)

    println("Starting analysis simulation for $t_max seconds")

    for i in 1:n_steps
        rk4_step!(dp, dt)

        times[i] = (i) * dt
        results[:, i] = dp.state
        history_theta1[i] = dp.state[1]
        history_theta2[i] = dp.state[2]

        energies[i] = total_energy(dp)

        # Store trajectory of the second mass
        x1, y1, x2, y2 = polar_to_cartesian(dp)
        traj_x[i] = x2
        traj_y[i] = y2
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

    data_to_save = [times energies energy_deviation]
    header = ["time_s" "energy_joules" "relative_error"]  # Columns: Time | Energy [J] | Relative Error

    open("./res/energy_stability.csv", "w") do io
        writedlm(io, header, ',')
        writedlm(io, data_to_save, ',')
    end

    p2 = plot(times, energy_deviation,
        label="Relative Error",
        xlabel="Time [s]", ylabel="ΔE / E0",
        title="Energy Stability (Precision)",
        color=:red, lw=1.5
    )

    savefig(p1, "./res/angles_evolution.png")
    savefig(p2, "./res/energy_stability.png")

    # Energy conservation statistics
    max_drift = maximum(abs.(energy_deviation))
    mean_drift = mean(abs.(energy_deviation))
    println("================================")
    println("Energy Conservation Statistics")
    println("Max drift: $(max_drift)")
    println("Mean drift: $(mean_drift)")
    println("================================")

    ### Generate GIF animation
    println("Generating GIF animation")

    # GIF parameters
    fps = 30                           # Frames per second for the animation
    step = Int(round(1 / (fps * dt)))  # Step between frames
    L = (dp.p1.l + dp.p2.l) * 1.1      # Graphical limits

    p = plot([0], [0], label="", color=:red, lw=1, alpha=0.5, # Trajectory in red
        xlims=(-L, L), ylims=(-L, L), aspect_ratio=:equal,
        grid=false, axis=false
    )
    plot!(p, [0, 0, 0], [0, 0, 0], color=:black, lw=3, marker=:circle, label="") # Pendulum

    anim = @animate for i in 1:step:n_steps
        print("\rRendering frame $i / $n_steps")

        # Update trajectory
        p[1][1][:x] = traj_x[1:i]
        p[1][1][:y] = traj_y[1:i]

        # Update pendulum position
        dp.state = results[:, i]
        x1, y1, x2, y2 = polar_to_cartesian(dp)
        p[1][2][:x] = [0, x1, x2]
        p[1][2][:y] = [0, y1, y2]

        # Update title
        title!(p, "Analysis (t = $(round(times[i], digits=1)) s)")
    end

    gif(anim, "./res/analysis_animation.gif", fps=fps)

    return times, results
end

function mode_comparison()
    """
    Compare the simulated motion of the double pendulum with experimental data.

    This function calculates initial conditions from experimental data, runs a simulation,
        and evaluates the root mean square error (RMSE) between the simulation and experimental results.
    """

    # Load data from CSV files
    data_path = "./res/video_data.csv"
    if !isfile(data_path)
        println("Error: data file not found.")
        return
    end

    # Read the CSV file
    data = readdlm(data_path, ';', skipstart=1)

    t_exp = data[:, 1]       # Time from experiment [s]
    theta1_exp = data[:, 2]  # Theta 1 from experiment [rad]
    theta2_exp = data[:, 3]  # Theta 2 from experiment [rad]

    # Unwrap experimental angles
    theta1_exp_unwrapped = zeros(length(theta1_exp))  # To store unwrapped theta 1 [rad]
    theta2_exp_unwrapped = zeros(length(theta2_exp))  # To store unwrapped theta 2 [rad]
    theta1_exp_unwrapped[1] = theta1_exp[1]           # Initialize first value
    theta2_exp_unwrapped[1] = theta2_exp[1]           # Initialize first value

    # Unwrapping process
    for i in 2:length(t_exp)
        theta1_exp_unwrapped[i] = unwrap_angle(theta1_exp[i], theta1_exp_unwrapped[i-1])
        theta2_exp_unwrapped[i] = unwrap_angle(theta2_exp[i], theta2_exp_unwrapped[i-1])
    end

    ### Initialize the simulation
    theta1_0 = theta1_exp_unwrapped[1]  # Initial angle theta 1 [rad]
    theta2_0 = theta2_exp_unwrapped[1]  # Initial angle theta 2 [rad]

    # Calculate initial angular velocities by averaging over N points
    N = min(10, length(t_exp) - 1)  # Use 10 points or fewer if not enough data

    # Angular velocity 1 (average of derivatives) [rad/s]
    omega1_0 = estimate_initial_angular_velocity(t_exp, theta1_exp_unwrapped, N)

    # Angular velocity 2 (average of derivatives) [rad/s]
    omega2_0 = estimate_initial_angular_velocity(t_exp, theta2_exp_unwrapped, N)

    initial_state = [theta1_0, theta2_0, omega1_0, omega2_0]
    dp = create_real_double_pendulum()
    dp.state = initial_state

    println("================================")
    println("Initial Conditions")
    display_dp_info(dp)
    println("================================")

    ### Simulation
    dt = 0.0001                       # Time step [s]
    t_max = maximum(t_exp)            # Simulation duration [s]
    n_steps = Int(floor(t_max / dt))  # Number of steps

    sim_time = zeros(n_steps)    # To store simulation time [s]
    sim_theta1 = zeros(n_steps)  # To store simulated theta 1 [rad]
    sim_theta2 = zeros(n_steps)  # To store simulated theta 2 [rad]

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
    title!(p1, "Theta 1: Reality vs Simulation")
    ylabel!(p1, "Angle [rad]")

    # Angle 2
    p2 = plot(t_exp, theta2_exp, label="Video", lw=2, color=:blue)
    plot!(p2, sim_time, wrap_angles(sim_theta2), label="Simulation", lw=2, color=:red)
    title!(p2, "Theta 2: Reality vs Simulation")
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

    # Calculate RMSE
    # Interpolation of simulation data at experimental time points
    sim_theta1_interp = [sim_theta1[argmin(abs.(sim_time .- t))] for t in t_exp]
    sim_theta2_interp = [sim_theta2[argmin(abs.(sim_time .- t))] for t in t_exp]

    println("================================")
    println("Errors (RMSE)")

    # Time segments [s]
    mask_p1 = t_exp .<= 1.0
    mask_p2 = (t_exp .> 1.0) .& (t_exp .<= 1.5)
    mask_p3 = t_exp .> 1.5

    function display_rmse_segment(mask, segment_name::String)
        """
        Calculate and display the RMSE for a specific time segment.

        Args:
            mask: A boolean mask indicating the time segment to analyze.
            segment_name: The name of the segment for display purposes.
        """
        if count(mask) == 0
            println("$segment_name: No data")
            return
        end

        # Local RMSE calculation
        err1 = sqrt(mean((sim_theta1_interp[mask] .- theta1_exp_unwrapped[mask]).^2))
        err2 = sqrt(mean((sim_theta2_interp[mask] .- theta2_exp_unwrapped[mask]).^2))

        # Formatted display
        println("$segment_name: θ1 = $(round(err1, digits=3)) rad, θ2 = $(round(err2, digits=3)) rad")
    end

    display_rmse_segment(mask_p1, "Phase 1 [0.0 - 1.0s]")
    display_rmse_segment(mask_p2, "Phase 2 [1.0 - 1.5s]")
    display_rmse_segment(mask_p3, "Phase 3 [> 1.5s]    ")

    # Total RMSE
    rmse_global1 = sqrt(mean((sim_theta1_interp .- theta1_exp_unwrapped).^2))
    rmse_global2 = sqrt(mean((sim_theta2_interp .- theta2_exp_unwrapped).^2))
    println("Global: θ1 = $(round(rmse_global1, digits=3)) rad, θ2 = $(round(rmse_global2, digits=3)) rad")
    println("================================")

    ### Generate comparative animation
    println("Generating comparative animation")

    dt_video = mean(diff(t_exp))    # Average time step from experimental data [s]
    fps_video = 1.0 / dt_video      # Frames per second in the video
    fps = 30                        # Desired frames per second for the animation
    speed_factor = fps / fps_video  # Speed adjustment factor

    L = (dp.p1.l + dp.p2.l) * 1.1   # Graphical limits

    # Plot initialisation
    p = plot([0], [0], label="", color=:blue, lw=1, alpha=0.3,
        xlims=(-L, L), ylims=(-L, L), aspect_ratio=:equal, grid=false, axis=false)                     # Real trail
    plot!(p, [0], [0], label="", color=:red, lw=1, alpha=0.3)                                          # Simulated trail
    plot!(p, [0, 0, 0], [0, 0, 0], label="Experimental", color=:blue, lw=3, alpha=0.6, marker=:circle) # Real pendulum
    plot!(p, [0, 0, 0], [0, 0, 0], label="Simulation", color=:red, lw=2, marker=:circle)               # Simulated pendulum

    # Store trails
    trail_real_x, trail_real_y = Float64[], Float64[]
    trail_sim_x, trail_sim_y = Float64[], Float64[]

    # Animation loop
    anim = @animate for i in 1:length(t_exp)
        print("\rRendering frame $i / $(length(t_exp))")

        # Update real pendulum
        dp.state[1] = theta1_exp_unwrapped[i]
        dp.state[2] = theta2_exp_unwrapped[i]
        x1_exp, y1_exp, x2_exp, y2_exp = polar_to_cartesian(dp)

        # Update real trail
        push!(trail_real_x, x2_exp)
        push!(trail_real_y, y2_exp)

        # Update simulated pendulum
        t_target = t_exp[i]
        idx_sim = argmin(abs.(sim_time .- t_target)) # Find closest simulation time index

        dp.state[1] = sim_theta1[idx_sim]
        dp.state[2] = sim_theta2[idx_sim]
        x1_sim, y1_sim, x2_sim, y2_sim = polar_to_cartesian(dp)

        # Update simulated trail
        push!(trail_sim_x, x2_sim)
        push!(trail_sim_y, y2_sim)

        # Update plot data
        # Update trails
        p[1][1][:x] = trail_real_x
        p[1][1][:y] = trail_real_y
        p[1][2][:x] = trail_sim_x
        p[1][2][:y] = trail_sim_y

        # Update pendulums
        p[1][3][:x] = [0, x1_exp, x2_exp]
        p[1][3][:y] = [0, y1_exp, y2_exp]
        p[1][4][:x] = [0, x1_sim, x2_sim]
        p[1][4][:y] = [0, y1_sim, y2_sim]

        title!(p, "Comparison (t = $(round(t_target, digits=1)) s) (x$(round(speed_factor, digits=2)))")
    end

    gif(anim, "./res/comparison_animation.gif", fps=fps)
end


function main()
    """
    Entry point of the program, allowing to choose between different operating modes.

    Modes:
        1. Real-time simulation.
        2. Analysis and results saving.
        3. Comparison with experimental data.
    """

    println("Choose the operating mode:")
    println("1. Real-time simulation")
    println("2. Analysis and results saving")
    println("3. Comparison with experimental data")
    print("Enter 1, 2 or 3: ")
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