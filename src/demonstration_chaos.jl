using Plots
using Statistics

include("utils.jl")

function calculate_lyapunov_exponent(times, diffs)
    """
    Calculate the Lyapunov exponent from divergences.
    Source: https://en.wikipedia.org/wiki/Lyapunov_exponent

    Args:
        times: Vector of time points [s]
        diffs: Vector of absolute divergences

    Returns:
        lambda: Lyapunov exponent [1/s]
    """

    # Convert to log(d(t)/d(0))
    log_diffs = log.(diffs ./ diffs[1])

    # Ignore the first 20% of points
    start_idx = Int(floor(length(times) * 0.2))
    t_fit = times[start_idx:end]
    log_fit = log_diffs[start_idx:end]

    # Linear regression: log(d/d0) = λ * t
    t_mean = mean(t_fit)
    log_mean = mean(log_fit)

    numerator = sum((t_fit .- t_mean) .* (log_fit .- log_mean))
    denominator = sum((t_fit .- t_mean).^2)

    lambda = numerator / denominator

    return lambda
end

function show_chaos_demo(epsilon::Float64=1e-6, duration::Float64=10.0)
    """
    Demonstration of Chaos in a Double Pendulum.
    Compare two nearly identical initial conditions and show their divergence over time.
    """

    dp1 = create_real_double_pendulum()
    dp2 = create_real_double_pendulum()
    dp2.state[2] += epsilon  # Small perturbation on theta 2

    dt = 0.0001                       # Time step [s]
    t_max = duration                  # Total duration [s]
    n_steps = Int(floor(t_max / dt))  # Number of steps

    times = zeros(n_steps)

    traj1_th1 = zeros(n_steps) # To store theta 1 of pendulum 1
    traj2_th1 = zeros(n_steps) # To store theta 1 of pendulum 2
    diffs_th1 = zeros(n_steps) # To store differences of theta 1

    traj1_th2 = zeros(n_steps) # To store theta 2 of pendulum 1
    traj2_th2 = zeros(n_steps) # To store theta 2 of pendulum 2
    diffs_th2 = zeros(n_steps) # To store differences of theta 2

    for i in 1:n_steps
        rk4_step!(dp1, dt)
        rk4_step!(dp2, dt)

        times[i] = i * dt

        # Store theta 1 and theta 2 trajectories and their differences
        traj1_th1[i] = dp1.state[1]
        traj2_th1[i] = dp2.state[1]
        diffs_th1[i] = abs(dp1.state[1] - dp2.state[1])

        traj1_th2[i] = dp1.state[2]
        traj2_th2[i] = dp2.state[2]
        diffs_th2[i] = abs(dp1.state[2] - dp2.state[2])
    end

    # Calculate Lyapunov exponents
    lambda_th1 = calculate_lyapunov_exponent(times, diffs_th1)
    lambda_th2 = calculate_lyapunov_exponent(times, diffs_th2)
    lambda_avg = (lambda_th1 + lambda_th2) / 2

    println("================================")
    println("Lyapunov Exponents")
    println("λ (theta1) = $(round(lambda_th1, digits=4)) s⁻¹")
    println("λ (theta2) = $(round(lambda_th2, digits=4)) s⁻¹")
    println("λ (average) = $(round(lambda_avg, digits=4)) s⁻¹")

    if lambda_avg > 0
        tau_lyapunov = 1.0 / lambda_avg
        t_double = log(2) / lambda_avg
        println("Lyapunov time: τ = $(round(tau_lyapunov, digits=2)) s")
        println("Doubling time: $(round(t_double, digits=2)) s")
    end
    println("================================")


    ### Plots
    # Theta 1 trajectories
    p1 = plot(times, [traj1_th1, traj2_th1],
        label=["Ref" "Perturbed"], color=[:blue :red], linestyle=[:solid :dash],
        title="Theta 1 Trajectories", ylabel="Angle [rad]", legend=:topleft)

    # Theta 2 trajectories
    p2 = plot(times, [traj1_th2, traj2_th2],
        label=["Ref" "Perturbed"], color=[:blue :red], linestyle=[:solid :dash],
        title="Theta 2 Trajectories", legend=:topleft)

    # Theta 1 Divergence (Log)
    p3 = plot(times, diffs_th1 .+ 1e-12, yscale=:log10,
        label="|Δθ1|", color=:purple,
        title="Divergence Log (Th1)", xlabel="Time [s]", ylabel="Log(|Δθ|)", legend=:topleft)

    # Add theoretical line e^(λt)
    theoretical_th1 = epsilon .* exp.(lambda_th1 .* times)
    plot!(p3, times, theoretical_th1,
        label="Theory: ε·e^(λt)", color=:orange, linestyle=:dash, lw=2)

    # Theta 2 Divergence (Log)
    p4 = plot(times, diffs_th2 .+ 1e-12, yscale=:log10,
        label="|Δθ2|", color=:purple,
        title="Divergence Log (Th2)", xlabel="Time [s]", legend=:topleft)

    # Add theoretical line
    theoretical_th2 = epsilon .* exp.(lambda_th2 .* times)
    plot!(p4, times, theoretical_th2,
        label="Theory: ε·e^(λt)", color=:orange, linestyle=:dash, lw=2)

    # Combine plots
    final_plot = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
    savefig(final_plot, "./res/demonstration_chaos.png")
end

show_chaos_demo()