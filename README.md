# Double Pendulum Simulation and Analysis

## Overview
This is a simulation and analysis project focused on the dynamics of a double pendulum. The project includes tools for real-time simulation, detailed analysis, and comparison with experimental data. It also demonstrates chaotic behavior in the double pendulum system.

## Features
- **Real-Time Simulation**: Visualize the motion of a double pendulum in real-time.
- **Analysis Mode**: Perform detailed simulations to study energy stability and angular evolution.
- **Comparison Mode**: Compare simulation results with experimental data.
- **Chaos Demonstration**: Illustrate the sensitivity to initial conditions using Lyapunov exponents.

## Project Structure
```
DoublePendulum/
├── res/                        
├── src/                        
│   ├── main.jl                 
│   ├── utils.jl                
│   ├── energies.jl             
│   ├── comparison.jl           
│   ├── demonstration_chaos.jl  
│   ├── system.jl               
│   ├── analytique.md           
│   └── extract_data.py         
├── test/                       
│   ├── test_energies.jl
│   ├── test_system.jl
│   └── test_utils.jl
└── README.md           
```

## Getting Started

### Prerequisites
- Julia (version 1.11.6 or higher)
- Python (for data extraction)
- Required Julia packages: `Plots`, `DelimitedFiles`, `Statistics`
- Required Python packages: `numpy`, `pandas`

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/louis272/DoublePendulum.git
   cd DoublePendulum
   ```
2. Install the required Julia packages:
   ```julia
   using Pkg
   Pkg.add(["Plots", "DelimitedFiles", "Statistics"])
   ```
3. Install the required Python packages:
   ```bash
   pip install numpy pandas
   ```

## Usage

### Running the Project
1. Open the `src/main.jl` file in your Julia environment.
2. Run the file to choose between the following modes:
   - **Real-Time Simulation**
   - **Analysis Mode**
   - **Comparison Mode**

### Chaos Demonstration
To demonstrate chaos, run the `src/demonstration_chaos.jl` file. This will generate plots illustrating the divergence of two nearly identical initial conditions.

### Optimization
To optimize the simulation parameters, you can use the `src/comparison.jl` script. 
It compares the simulation results with experimental data and adjusts parameters to minimize discrepancies.

### Outputs
- Simulation results and plots are saved in the `res/` directory.
- Example outputs include:
  - `angles_evolution.png`: Angular evolution over time.
  - `energy_stability.png`: Energy stability analysis.
  - `comparison_reality_simulation_wrapped.png`: Comparison of wrapped angles.
  - `demonstration_chaos.png`: Chaos demonstration plots.

## Analytical Part
For the analytical part of the project, refer to the `src/analytique.md` file.
It contains detailed derivations and explanations of the equations governing the double pendulum system.