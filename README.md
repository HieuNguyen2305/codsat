# Energy Efficiency Optimization in Integrated Satellite-Terrestrial Networks

This repository implements optimization algorithms to maximize downlink Energy Efficiency (EE) in Integrated Satellite-Terrestrial Networks (ISTN). The project focuses on solving non-convex fractional programming problems using advanced optimization techniques.

## 🚀 Key Features & Techniques
The solution leverages the following mathematical methods:
* **Dinkelbach's Algorithm:** To transform the fractional objective function (Energy Efficiency) into a subtractive parametric form.
* **SCA (Successive Convex Approximation):** To handle non-convex constraints by approximating them with convex functions.
* **ADMM (Alternating Direction Method of Multipliers):** Implemented for distributed optimization/resource allocation.

## 📂 File Structure

### Main Scripts
* `main.m`: The primary entry point for the centralized optimization scheme.
* `admm_main.m`: The main script for running the ADMM-based distributed algorithm.
* `fivescheme.m`: Implementation of 5 different baseline schemes for comparison.

### Optimization Modules
* `dinkelbach_sca.m`: Core function combining Dinkelbach and SCA.
* `admm_dinkelbach.m`: ADMM implementation integrated with Dinkelbach's method.
* `sca_subproblem.m`: Solver for the SCA sub-problems.
* `solve_bs_local.m`: Local optimization solver for the Base Station (BS).
* `solve_st_local.m`: Local optimization solver for the Satellite (ST).

### Utilities & Data
* `params.m`: System configuration and simulation parameters.
* `gen_channels.m`: Channel generation and modeling.
* `find_init_point.m`: Helper function to find feasible initialization points.
* `sosanh.m`: Script for plotting and comparing results.
* `compare.xlsx`: Excel file storing simulation results for analysis.

## 🛠 Requirements
* **MATLAB** (R2020b or later recommended)
* **CVX** (Matlab Software for Disciplined Convex Programming) - *Assuming you are using CVX based on the optimization context.*

## usage
1.  Configure the simulation parameters in `params.m`.
2.  Run `main.m` to execute the proposed algorithm.
3.  Run `admm_main.m` to test the distributed approach.
4.  Use `sosanh.m` to visualize the performance comparison.

## 👨‍💻 Author
**Hieu Nguyen** (HieuNguyen2305)
