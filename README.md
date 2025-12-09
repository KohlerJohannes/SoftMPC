# **A model predictive control framework with robust stability under arbitrarily large disturbances** - Simulation Code

This repository contains the MATLAB code that accompanies the paper

Johannes Köhler,  Melanie Zeilinger "A model predictive control framework with robust stability guarantees under unbounded disturbances", [ArXiv link](https://arxiv.org/abs/2207.10216).

# Dependencies
- [YALMIP](https://yalmip.github.io/): formulate linear MPC problems
-  [MPT-3 toolbox](https://www.mpt3.org/): plot polytopic sets
- [CasADi](https://web.casadi.org/) + [IPOPT](https://coin-or.github.io/Ipopt/) solve nonlinear MPC problems


# Linear_MassSrping
 - **main.m** runs the offline design, closed-loop simulations and computational timings

 - **analyse_results.m** generates plots and quantitative results (comp. times, performance, etc.)

# Nonlinear_FourTank
 
 - **main.m** runs the offline design, closed-loop simulations and computational timings
             

 - **analyse_results.m** generates plots and lists comp. times


 
