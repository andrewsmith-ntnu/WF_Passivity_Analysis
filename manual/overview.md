# Folder overview

This page lists the functions and scripts in the Public folder with a one-line summary. For full per-file details see `functions.md` and `scripts.md`.

## Functions

- **IM_C.m** — Sequence-domain impedance for a capacitor
- **IM_RL.m** — Sequence-domain impedance of an RL element
- **IM_WT.m** — Compute sequence-domain impedance FRD for a WT converter model
- **NR_steadystate_VSC_GF.m** — Calculate the steady state operating point of a grid-following VSC using Newton-Raphson method.
- **base_values.m** — Compute per-unit base values for a converter struct
- **calculate_WF.m** — Calculates the wind farm parameters and impedance models based on the wind turbine and transmission system data.
- **function_NR_VSC_GF.m** — Calculates the Newton-Raphson equations and Jacobian for the VSC grid-forming converter model.
- **impedance_WT_agg.m** — Aggregate single-WT model into equivalent aggregated model
- **impedance_WT_blackbox.m** — Build black-box WT aggregated impedance model
- **impedance_WT_blackbox_obf.m** — Calculates the state-space model and impedance of a black-boxed wind turbine converter
- **initialize_NR_VSC_GF.m** — Initializes the state variables and steady-state values for the Newton-Raphson VSC grid-forming converter model.
- **initialize_WT_agg.m** — Initializes the aggregated wind turbine converter model based on the single wind turbine converter parameters and aggregated grid-side impedance.
- **initialize_WT_default_gfl.m** — Initializes the wind turbine converter model in grid-following mode with default parameters.
- **initialize_WT_default_gfm.m** — Initializes the wind turbine converter model in grid-forming mode with default parameters.
- **matrix_VSC_GF.m** — Calculate the matrixes for the linear model of the VSC grid-following converter.
- **matrix_VSM_DEM.m** — Calculate the matrixes for the linear model of the VSM DEM converter with PLL.
- **matrix_VSM_DEM_noPLL.m** — Calculate the matrixes for the linear model of the VSM DEM converter without PLL.
- **matrix_VSM_QSEM.m** — Calculate the matrixes for the linear model of the VSM QSEM converter with PLL.
- **matrix_VSM_QSEM_noPLL.m** — Calculate the matrixes for the linear model of the VSM QSEM converter without PLL.
- **matrix_VSM_VC.m** — Calculate the matrixes for the linear model of the VSM VC converter with PLL.
- **matrix_VSM_VC_noPLL.m** — Calculate the matrixes for the linear model of the VSM VC converter without PLL.
- **mimo_passivity.m** — Simple positive-definiteness check for 2x2 FRD impedance
- **plot_p_vary.m** — Plot passivity and impedance variation across datasets
- **plot_pass_imp.m** — Plot passivity and impedance magnitude/imaginary part
- **setNestedField.m** — Set nested struct field using dot-separated name
- **solveNewtonRaphson.m** — Solve a system of nonlinear equations using the Newton-Raphson method with backtracking.
- **steadystate_VSC_GF.m** — Calculate the steady state operating point of a grid-following VSC.
- **steadystate_VSM_DEM.m** — Calculate the steady state operating point of a grid-forming VSM DEM converter with PLL
- **steadystate_VSM_DEM_noPLL.m** — Calculate the steady state operating point of a grid-forming VSM DEM converter without PLL
- **steadystate_VSM_DEM_no_PLL.m** — Calculate the steady state operating point of a grid-forming VSM DEM converter without PLL
- **steadystate_VSM_QSEM.m** — Calculate the steady state operating point of a grid-forming VSM QSEM converter with PLL
- **steadystate_VSM_QSEM_noPLL.m** — Calculate the steady state operating point of a grid-forming VSM QSEM converter without PLL
- **steadystate_VSM_VC.m** — Calculate the steady state operating point of a grid-forming VSM VC converter with PLL
- **steadystate_VSM_VC_noPLL.m** — Calculate the steady state operating point of a grid-forming VSM VC converter without PLL
- **steadystate_nonlinear_VSC_GF.m** — Steady-state nonlinear equations for a grid-following VSC
- **steadystate_nonlinear_VSM_DEM.m** — Steady-state nonlinear equations for a grid-forming VSM DEM converter with PLL
- **steadystate_nonlinear_VSM_DEM_noPLL.m** — Steady-state nonlinear equations for a grid-forming VSM DEM converter without PLL
- **steadystate_nonlinear_VSM_QSEM.m** — Steady-state nonlinear equations for a grid-forming VSM QSEM converter with PLL
- **steadystate_nonlinear_VSM_QSEM_noPLL.m** — Steady-state nonlinear equations for a grid-forming VSM QSEM converter without PLL
- **steadystate_nonlinear_VSM_VC.m** — Steady-state nonlinear equations for a grid-forming VSM VC converter with PLL
- **steadystate_nonlinear_VSM_VC_noPLL.m** — Steady-state nonlinear equations for a grid-forming VSM VC converter without PLL

## Scripts

- **WF_parameter_variation.m** — Main script to run parameter variation cases for WF impedance modeling.
- **parameters_WF_default.m** — This script initializes default parameters for wind farm (WF) modeling. Use this to change the default values.

