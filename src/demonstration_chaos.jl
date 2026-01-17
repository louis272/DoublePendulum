using Plots
using Statistics

include("utils.jl")

function show_chaos_demo(epsilon::Float64=1e-6, duration::Float64=10.0)
    """
    Demonstration of Chaos in a Double Pendulum.
    Compares two nearly identical initial conditions and shows their divergence over time.
    """

    dp1 = create_real_double_pendulum()
    dp2 = create_real_double_pendulum()
    dp2.state[2] += epsilon  # Small perturbation on theta2

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
        title="Divergence Log (Th1)", xlabel="Temps [s]", ylabel="Log(|Δθ|)")

    # Theta 2 Divergence (Log)
    p4 = plot(times, diffs_th2 .+ 1e-12, yscale=:log10,
        label="|Δθ2|", color=:purple,
        title="Divergence Log (Th2)", xlabel="Temps [s]")

    # Combine plots
    final_plot = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
    savefig(final_plot, "./res/demonstration_chaos.png")
end

show_chaos_demo()