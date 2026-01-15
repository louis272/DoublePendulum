using Plots
using DelimitedFiles
using Statistics

include("utils.jl")

function compare_lengths(l1_measure_mm::Float64, l2_measure_mm::Float64, csv_path::String)
    """
    Compares the lengths of the bars measured by hand (in mm) with those extracted from the video (in pixels).
    This function allows to:
    1. Calculate the scale (pixels/mm) of the video.
    2. Check the stability of the tracking (the noise on the length).
    3. Detect perspective problems (if the scale calculated for L1 is different from L2).

    Args :
        l1_measure_mm : Measured length of the first bar (mm).
        l2_measure_mm : Measured length of the second bar (mm).
        csv_path : Path to the CSV file containing the video data.
    """

    # Load data
    if !isfile(csv_path)
        println("Error: File not found : $csv_path")
        return
    end

    # Read CSV file
    data = readdlm(csv_path, ';', skipstart=1)

    time = data[:, 1]
    l1_px = data[:, 4] # Length L1 in pixels
    l2_px = data[:, 5] # Length L2 in pixels

    # Statistical Analysis
    mean_l1_px = mean(l1_px)
    mean_l2_px = mean(l2_px)
    std_l1 = std(l1_px)
    std_l2 = std(l2_px)

    # Calculation of the scale (Calibration) using L1 as a reference (since attached to the fixed pivot)
    scale_px_per_mm = mean_l1_px / l1_measure_mm

    # Convert video data to mm using this scale
    l1_video_mm = l1_px ./ scale_px_per_mm
    l2_video_mm = l2_px ./ scale_px_per_mm

    # Comparison of scales (Perspective test)
    scale_l2 = mean_l2_px / l2_measure_mm
    ecart_echelle = abs(scale_px_per_mm - scale_l2) / scale_px_per_mm * 100

    println("LENGTH COMPARISON REPORT")
    println("Estimated scale (via L1) : $(round(scale_px_per_mm, digits=2)) px/mm")
    println("Tracking stability L1 : ± $(round(std_l1, digits=2)) px")
    println("Tracking stability L2 : ± $(round(std_l2, digits=2)) px")
    println("Consistency L1 vs L2 : Scale difference of $(round(ecart_echelle, digits=2)) %")

    # Graphs
    # Graph L1
    p1 = plot(time, l1_video_mm, label="L1 Video (reconstructed)", lw=2, color=:blue)
    hline!(p1, [l1_measure_mm], label="L1 Measured ($l1_measure_mm mm)", color=:black, linestyle=:dash, lw=2)
    ylabel!(p1, "Length [mm]")
    title!(p1, "Stability L1 (Pivot -> Mass 1)")
    ylims!(p1, l1_measure_mm * 0.9, l1_measure_mm * 1.1) # Zoom ±10%

    # Graph L2
    p2 = plot(time, l2_video_mm, label="L2 Video (reconstructed)", lw=2, color=:red)
    hline!(p2, [l2_measure_mm], label="L2 Measured ($l2_measure_mm mm)", color=:black, linestyle=:dash, lw=2)
    xlabel!(p2, "Time [s]")
    ylabel!(p2, "Length [mm]")
    title!(p2, "Stability L2 (Mass 1 -> Mass 2)")
    ylims!(p2, l2_measure_mm * 0.9, l2_measure_mm * 1.1)

    final_plot = plot(p1, p2, layout=(2, 1), size=(800, 600))
    savefig(final_plot, "./res/lengths_comparison_$l2_measure_mm.png")
end

function optimize_l2(l1_ref_mm::Float64, l2_guess_mm::Float64, csv_path::String; range_percent=10.0, steps=200)
    """
    Finds the value of L2 that minimizes the scale inconsistency between the two bars.

    Args :
        l1_ref_mm : The reliable reference length (Bar 1 measured by hand) [mm].
        l2_guess_mm : An estimate of bar 2 (to center the search) [mm].
        range_percent : Search range around the estimate (by default ±10%).
        steps : Number of steps in the scanning method (by default 200).

    Returns :
        The optimal value of L2 [mm] that minimizes the scale inconsistency.
    """

    # Load and read data
    if !isfile(csv_path)
        println("Error: File not found : $csv_path")
        return
    end
    data = readdlm(csv_path, ';', skipstart=1)

    l1_px_mean = mean(data[:, 4])
    l2_px_mean = mean(data[:, 5])

    # Direct Calculation
    # If Scale = l1_px / l1_mm, then l2_mm = l2_px / Scale
    scale_ref = l1_px_mean / l1_ref_mm
    l2_optimal_direct = l2_px_mean / scale_ref

    println("L2 OPTIMIZATION")
    println("Mean L1 (pixels) : $(round(l1_px_mean, digits=2))")
    println("Mean L2 (pixels) : $(round(l2_px_mean, digits=2))")
    println("Detected scale    : $(round(scale_ref, digits=2)) px/mm")
    println("------------------------------------------------")
    println("Direct suggestion  : L2 = $(round(l2_optimal_direct, digits=4)) mm")

    # Scanning Method (To visualize the error dip)
    l2_min = l2_guess_mm * (1 - range_percent/100)
    l2_max = l2_guess_mm * (1 + range_percent/100)
    l2_values = range(l2_min, l2_max, length=steps)
    errors = zeros(steps)

    for (i, l2_test) in enumerate(l2_values)
        # Calculate the scale deviation for this candidate l2_test
        scale_test = l2_px_mean / l2_test
        # Relative error between L1 scale and L2 scale (%)
        errors[i] = abs(scale_ref - scale_test) / scale_ref * 100
    end

    # Find the minimum of the scan
    min_err, idx_min = findmin(errors)
    l2_best_scan = l2_values[idx_min]

    # Convergence graph
    p = plot(l2_values, errors,
        label="Consistency Error (%)",
        xlabel="Tested L2 Value (mm)",
        ylabel="Scale Deviation (%)",
        title="L2 Length Optimization",
        lw=2, color=:blue
    )
    vline!(p, [l2_best_scan], label="Detected Optimum ($(round(l2_best_scan, digits=2)) mm)", color=:red, linestyle=:dash)

    savefig(p, "./res/optimization_l2_$l2_guess_mm.png")

    return l2_optimal_direct
end

function simulate_double_pendulum(m1, m2, times)
    dp = create_real_double_pendulum()

    dp.p1.m = m1
    dp.p2.m = m2

    t_max = 2.0 # Simulation duration [s]
    n_steps = length(times)

    ### Store results
    results = zeros(4, n_steps)

    for i in 1:n_steps
        dt = times[i] - (i > 1 ? times[i-1] : 0.0) # Time step [s]
        rk4_step!(dp, dt)
        results[:, i] = dp.state  # [theta1_0, theta2_0, omega1_0, omega2_0] -> one column per time step
    end

    return results
end

function calculate_error(simulated_data, experimental_data)
    """
    Calculate the error between simulated and experimental data.

    Args :
        simulated_data : The simulated data from the model.
        experimental_data : The experimental data from the CSV file.

    Returns :
        A scalar error value representing the discrepancy.
    """

    # Use the sum of squared differences for theta1 and theta2
    sim_theta1 = simulated_data[1, :]
    sim_theta2 = simulated_data[2, :]

    exp_theta1 = experimental_data[:, 2]
    exp_theta2 = experimental_data[:, 3]

    error_theta1 = sum((sim_theta1 .- exp_theta1).^2)
    error_theta2 = sum((sim_theta2 .- exp_theta2).^2)

    return error_theta1, error_theta2
end

function optimize_m1_m2(m1_ref::Float64, m2_ref::Float64, csv_path::String, range_percent=10.0, steps=200)
    """
    Finds the optimal values of masses m1 and m2 that minimize the discrepancy between simulated and experimental data.

    Args :
        m1_ref : Reference mass 1 [kg].
        m2_ref : Reference mass 2 [kg].
        csv_path : Path to the CSV file containing the video data.
        range_percent : Search range around the reference masses (by default ±10%).
        steps : Number of steps in the scanning method (by default 200).

    Returns :
        Two tuples (m1_optimal, m2_optimal) representing the optimal masses [kg] for theta1 and theta2.
    """

    # Load and read data
    if !isfile(csv_path)
        println("Error: File not found : $csv_path")
        return
    end
    data = readdlm(csv_path, ';', skipstart=1)

    # Perform optimization by simulating the double pendulum with various mass combinations and comparing to the experimental data.
    errors_theta1 = zeros(steps, steps)
    errors_theta2 = zeros(steps, steps)
    m1_values = range(m1_ref * (1 - range_percent / 100), m1_ref * (1 + range_percent / 100), length=steps)
    m2_values = range(m2_ref * (1 - range_percent / 100), m2_ref * (1 + range_percent / 100), length=steps)

    for (i, m1_test) in enumerate(m1_values)
        for (j, m2_test) in enumerate(m2_values)
            # Simulate the system with the current mass combination
            simulated_data = simulate_double_pendulum(m1_test, m2_test, data[:, 1])
            # Calculate the error between simulated and experimental data
            errors_theta1[i, j], errors_theta2[i, j] = calculate_error(simulated_data, data)
        end
    end

    # Find the combination of masses that minimizes the error for theta 1
    min_err, idx = findmin(errors_theta1)
    m1_optimal_theta1 = m1_values[idx[1]]
    m2_optimal_theta1 = m2_values[idx[2]]

    # Find the combination of masses that minimizes the error for theta 2
    min_err, idx = findmin(errors_theta2)
    m1_optimal_theta2 = m1_values[idx[1]]
    m2_optimal_theta2 = m2_values[idx[2]]

    return (m1_optimal_theta1, m2_optimal_theta1), (m1_optimal_theta2, m2_optimal_theta2)
end


# compare_lengths(91.74, 69.33, "./res/video_data.csv")
#
# l2_optimal = optimize_l2(91.74, 69.33, "./res/video_data.csv")
# println("Value retained for simulation : L2 = $(round(l2_optimal, digits=2))")
#
# compare_lengths(91.74, round(l2_optimal, digits=2), "./res/video_data.csv")

(m1_optimal_theta1, m2_optimal_theta1), (m1_optimal_theta2, m2_optimal_theta2) = optimize_m1_m2(30.00e-3, 2.00e-3, "./res/video_data.csv")
println("Optimal masses for θ1 fit : m1 = $(round(m1_optimal_theta1*1e3, digits=3)) g, m2 = $(round(m2_optimal_theta1*1e3, digits=3)) g")
println("Optimal masses for θ2 fit : m1 = $(round(m1_optimal_theta2*1e3, digits=3)) g, m2 = $(round(m2_optimal_theta2*1e3, digits=3)) g")