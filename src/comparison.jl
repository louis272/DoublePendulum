using Plots
using DelimitedFiles
using Statistics

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

compare_lengths(91.74, 69.33, "./res/video_data.csv")

l2_optimal = optimize_l2(91.74, 69.33, "./res/video_data.csv")
println("Value retained for simulation : L2 = $(round(l2_optimal, digits=2))")

compare_lengths(91.74, round(l2_optimal, digits=2), "./res/video_data.csv")
