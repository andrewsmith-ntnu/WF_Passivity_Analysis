# Functions in C:\Users\andrsmit\OneDrive - NTNU\Papers\WF impedance aggregation\Public

### IM_C.m

**Signature:** `function [Z_frd]=IM_C(f,f1,C)`

**Summary:** Sequence-domain impedance for a capacitor

#### Help

```matlab
Sequence-domain impedance for a capacitor

[Z_frd] = IM_C(f, f1, C)

Inputs
  f  - frequency vector (Hz)
  f1 - fundamental frequency (Hz)
  C  - capacitance (F)

Outputs
  Z_frd - frequency response data object (2x2 sequence-domain impedance)

```

### IM_RL.m

**Signature:** `function [Z_frd]=IM_RL(f,f1,R,L)`

**Summary:** Sequence-domain impedance of an RL element

#### Help

```matlab
Sequence-domain impedance of an RL element

[Z_frd] = IM_RL(f, f1, R, L)

Inputs
  f  - frequency vector (Hz)
  f1 - fundamental frequency (Hz)
  R  - resistance (Ohm)
  L  - inductance (H)

Outputs
  Z_frd - FRD object containing 2x2 sequence-domain impedance (pos/neg)

```

### IM_WT.m

**Signature:** `function [Z_frd]=IM_WT(f,cvtr)`

**Summary:** Compute sequence-domain impedance FRD for a WT converter model

#### Help

```matlab
Compute sequence-domain impedance FRD for a WT converter model

[Z_frd] = IM_WT(f, cvtr)

Inputs
  f    - frequency vector (Hz)
  cvtr - converter struct containing state-space matrices A,B,C and parameters

Outputs
  Z_frd - frequency response data (frd) object of the sequence-domain impedance

Notes
  This function builds a state-space model from cvtr.A/B/C, selects the
  relevant voltage/current channels, computes the impedance in dq and
  transforms it to positive/negative sequence domain.


Build continuous-time state-space system from provided matrices
```

### NR_steadystate_VSC_GF.m

**Signature:** `function [cvtr] = NR_steadystate_VSC_GF(cvtr)`

**Summary:** Calculate the steady state operating point of a grid-following VSC using Newton-Raphson method.

#### Help

```matlab
Calculate the steady state operating point of a grid-following VSC using Newton-Raphson method.

% 'solveNewtonRaphson' function parameters
```

### base_values.m

**Signature:** `function [cvtr] = base_values(cvtr)`

**Summary:** Compute per-unit base values for a converter struct

#### Help

```matlab
Compute per-unit base values for a converter struct

cvtr = base_values(cvtr)

Inputs
  cvtr - struct with fields 'parameters' containing Vnr, Inr, fnr, etc.

Outputs
  cvtr - same struct with populated cvtr.pu.* fields used throughout the code

Base apparent power (kVA)
```

### calculate_WF.m

**Signature:** `function WF = calculate_WF(WF)`

**Summary:** Calculates the wind farm parameters and impedance models based on the wind turbine and transmission system data.

#### Help

```matlab
Calculates the wind farm parameters and impedance models based on the wind turbine and transmission system data.

```

### function_NR_VSC_GF.m

**Signature:** `function[y,dy]=function_NR_VSC_GF(x,conv)`

**Summary:** Calculates the Newton-Raphson equations and Jacobian for the VSC grid-forming converter model.

#### Help

```matlab
Calculates the Newton-Raphson equations and Jacobian for the VSC grid-forming converter model.

```

### impedance_WT_agg.m

**Signature:** `function [cvtr_agg,Z_cvtr_agg] = impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg)`

**Summary:** Aggregate single-WT model into equivalent aggregated model

#### Help

```matlab
Aggregate single-WT model into equivalent aggregated model

[cvtr_agg, Z_cvtr_agg] = impedance_WT_agg(freq, cvtr, Sb_WF, r_g_agg, l_g_agg)

Inputs
  freq    - frequency vector (Hz) for impedance calculation
  cvtr    - base converter model struct
  Sb_WF   - apparent power of the wind farm (per unit or absolute depending on usage)
  r_g_agg - equivalent grid resistance for aggregated model
  l_g_agg - equivalent grid inductance for aggregated model

Outputs
  cvtr_agg    - converter struct updated for aggregated ratings and parameters
  Z_cvtr_agg  - aggregated converter impedance (sequence-domain FRD)


```

### impedance_WT_blackbox.m

**Signature:** `function [cvtr_ss,Z_cvtr_agg] = impedance_WT_blackbox(freq,f1,V_wt,V_coll,Sb_WT,n_t,ref,rg_eq,lg_eq)`

**Summary:** Build black-box WT aggregated impedance model

#### Help

```matlab
Build black-box WT aggregated impedance model

[cvtr_ss, Z_cvtr_agg] = impedance_WT_blackbox(freq, f1, V_wt, V_coll, Sb_WT, n_t, ref, rg_eq, lg_eq)

Inputs
  freq   - frequency vector (Hz)
  f1     - fundamental frequency (Hz)
  V_wt   - nominal winding voltage of a WT (V)
  V_coll - collection voltage used for scaling
  Sb_WT  - rated apparent power of one WT
  n_t    - number of turbines in the farm
  ref    - references struct with fields pref, qref, vref, omegaref
  rg_eq  - equivalent grid resistance for black-box aggregation
  lg_eq  - equivalent grid inductance for black-box aggregation

Outputs
  cvtr_ss     - steady-state converter struct used for black-box model
  Z_cvtr_agg  - aggregated impedance (sequence-domain FRD)

```

### impedance_WT_blackbox_obf.m

**Signature:** `function [cvtr_ss,Z_cvtr_agg] = impedance_WT_blackbox_obf(freq,f1,V_wt,V_coll,Sb_WT,n_t,ref,rg_eq,lg_eq)`

**Summary:** Calculates the state-space model and impedance of a black-boxed wind turbine converter

#### Help

```matlab
Calculates the state-space model and impedance of a black-boxed wind turbine converter
```

### initialize_NR_VSC_GF.m

**Signature:** `function conv = initialize_NR_VSC_GF(conv)`

**Summary:** Initializes the state variables and steady-state values for the Newton-Raphson VSC grid-forming converter model.

#### Help

```matlab
Initializes the state variables and steady-state values for the Newton-Raphson VSC grid-forming converter model.

```

### initialize_WT_agg.m

**Signature:** `function [cvtr_agg] = initialize_WT_agg(cvtr,Sb_WF,r_g_agg,l_g_agg)`

**Summary:** Initializes the aggregated wind turbine converter model based on the single wind turbine converter parameters and aggregated grid-side impedance.

#### Help

```matlab
Initializes the aggregated wind turbine converter model based on the single wind turbine converter parameters and aggregated grid-side impedance.

```

### initialize_WT_default_gfl.m

**Signature:** `function cvtr = initialize_WT_default_gfl(WF)`

**Summary:** Initializes the wind turbine converter model in grid-following mode with default parameters.

#### Help

```matlab
Initializes the wind turbine converter model in grid-following mode with default parameters.

CC_DEM_noPLL = 0, CC_DEM = 1, CC_QSEM_noPLL = 2, CC_QSEM =3, VCVSM_noPLL, = 4, VCVSM = 5, GF-PQ = 6

```

### initialize_WT_default_gfm.m

**Signature:** `function cvtr = initialize_WT_default_gfm(WF)`

**Summary:** Initializes the wind turbine converter model in grid-forming mode with default parameters.

#### Help

```matlab
Initializes the wind turbine converter model in grid-forming mode with default parameters.

CC_DEM_noPLL = 0, CC_DEM = 1, CC_QSEM_noPLL = 2, CC_QSEM =3, VCVSM_noPLL, = 4, VCVSM = 5, GFL-PQ = 6

```

### matrix_VSC_GF.m

**Signature:** `function conv = matrix_VSC_GF(conv,conv_str)`

**Summary:** Calculate the matrixes for the linear model of the VSC grid-following converter.

#### Help

```matlab
Calculate the matrixes for the linear model of the VSC grid-following converter.

% Calculate the matrixes for the linear model
```

### matrix_VSM_DEM.m

**Signature:** `function conv = matrix_VSM_DEM(conv,conv_str)`

**Summary:** Calculate the matrixes for the linear model of the VSM DEM converter with PLL.

#### Help

```matlab
Calculate the matrixes for the linear model of the VSM DEM converter with PLL.

% Calculate the matrixes for the linear model
```

### matrix_VSM_DEM_noPLL.m

**Signature:** `function conv = matrix_VSM_DEM_noPLL(conv,conv_str)`

**Summary:** Calculate the matrixes for the linear model of the VSM DEM converter without PLL.

#### Help

```matlab
Calculate the matrixes for the linear model of the VSM DEM converter without PLL.

% Calculate the matrixes for the linear model
```

### matrix_VSM_QSEM.m

**Signature:** `function conv = matrix_VSM_QSEM(conv,conv_str)`

**Summary:** Calculate the matrixes for the linear model of the VSM QSEM converter with PLL.

#### Help

```matlab
Calculate the matrixes for the linear model of the VSM QSEM converter with PLL.

% Calculate the matrixes for the linear model
```

### matrix_VSM_QSEM_noPLL.m

**Signature:** `function conv = matrix_VSM_QSEM_noPLL(conv,conv_str)`

**Summary:** Calculate the matrixes for the linear model of the VSM QSEM converter without PLL.

#### Help

```matlab
Calculate the matrixes for the linear model of the VSM QSEM converter without PLL.

% Calculate the matrixes for the linear model
```

### matrix_VSM_VC.m

**Signature:** `function conv = matrix_VSM_VC(conv,conv_str)`

**Summary:** Calculate the matrixes for the linear model of the VSM VC converter with PLL.

#### Help

```matlab
Calculate the matrixes for the linear model of the VSM VC converter with PLL.

% Calculate the matrixes for the linear model
```

### matrix_VSM_VC_noPLL.m

**Signature:** `function conv = matrix_VSM_VC_noPLL(conv,conv_str)`

**Summary:** Calculate the matrixes for the linear model of the VSM VC converter without PLL.

#### Help

```matlab
Calculate the matrixes for the linear model of the VSM VC converter without PLL.

% Calculate the matrixes for the linear model
```

### mimo_passivity.m

**Signature:** `function pass = mimo_passivity(Z)`

**Summary:** Simple positive-definiteness check for 2x2 FRD impedance

#### Help

```matlab
Simple positive-definiteness check for 2x2 FRD impedance

pass = mimo_passivity(Z)

Inputs
  Z - 2x2 FRD object representing impedance (Z.ResponseData is 2x2xN)

Outputs
  pass - vector containing the minimum eigenvalue of the Hermitian part
         (used as a passivity indicator) for each frequency sample

Notes
  This routine forms the Hermitian part P = Z + Z^H and computes eigenvalues
  per frequency. Positive minimum eigenvalue suggests passivity.

Extract 2x2 FRD blocks
```

### plot_p_vary.m

**Signature:** `function [] = plot_p_vary(f, Z1_array, Z2_array, n_o, n_i, opts)`

**Summary:** Plot passivity and impedance variation across datasets

#### Help

```matlab
Plot passivity and impedance variation across datasets

plot_p_vary(f, Z1_array, Z2_array, n_o, n_i, opts)

Inputs
  f         - figure handle or figure index
  Z1_array  - cell array of FRD objects used to evaluate passivity
  Z2_array  - cell array of FRD objects to plot impedance magnitude/imag
  n_o, n_i  - output/input indices (for selecting matrix element)
  opts      - optional struct: fields 'plot_pass' (bool), 'plot_type' ('imag')

This utility plots passivity (min eigenvalue of Hermitian part) and
either |Z| or Im(Z) for several parameter sweeps contained in Z2_array.

```

### plot_pass_imp.m

**Signature:** `function [] = plot_pass_imp(f, Z1, Z2, n_o, n_i, plot_type)`

**Summary:** Plot passivity and impedance magnitude/imaginary part

#### Help

```matlab
Plot passivity and impedance magnitude/imaginary part

plot_pass_imp(f, Z1, Z2, n_o, n_i, plot_type)

Inputs
  f         - figure handle or index
  Z1        - FRD used for passivity calculation
  Z2        - FRD whose element is plotted
  n_o, n_i  - output/input indices for selecting element from Z2
  plot_type - optional, 'imag' to plot imaginary part instead of magnitude

```

### setNestedField.m

**Signature:** `function S = setNestedField(S, fieldName, value)`

**Summary:** Set nested struct field using dot-separated name

#### Help

```matlab
Set nested struct field using dot-separated name

S = setNestedField(S, 'a.b.c', value)

This utility uses subsasgn to assign a nested field specified as a
dot-separated string. It creates the appropriate substruct indexing.

```

### solveNewtonRaphson.m

**Signature:** `function [xs,dxs] = solveNewtonRaphson(fun,ys,xg,cnv,Ns,bta,litmax)`

**Summary:** Solve a system of nonlinear equations using the Newton-Raphson method with backtracking.

#### Help

```matlab
Solve a system of nonlinear equations using the Newton-Raphson method with backtracking.

```

### steadystate_VSC_GF.m

**Signature:** `function conv = steadystate_VSC_GF(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-following VSC.

#### Help

```matlab
Calculate the steady state operating point of a grid-following VSC.

```

### steadystate_VSM_DEM.m

**Signature:** `function conv = steadystate_VSM_DEM(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-forming VSM DEM converter with PLL

#### Help

```matlab
Calculate the steady state operating point of a grid-forming VSM DEM converter with PLL

```

### steadystate_VSM_DEM_noPLL.m

**Signature:** `function conv = steadystate_VSM_DEM_noPLL(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-forming VSM DEM converter without PLL

#### Help

```matlab
Calculate the steady state operating point of a grid-forming VSM DEM converter without PLL

```

### steadystate_VSM_DEM_no_PLL.m

**Signature:** `function conv = steadystate_VSM_DEM_no_PLL(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-forming VSM DEM converter without PLL

#### Help

```matlab
Calculate the steady state operating point of a grid-forming VSM DEM converter without PLL

```

### steadystate_VSM_QSEM.m

**Signature:** `function conv = steadystate_VSM_QSEM(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-forming VSM QSEM converter with PLL

#### Help

```matlab
Calculate the steady state operating point of a grid-forming VSM QSEM converter with PLL

```

### steadystate_VSM_QSEM_noPLL.m

**Signature:** `function conv = steadystate_VSM_QSEM_noPLL(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-forming VSM QSEM converter without PLL

#### Help

```matlab
Calculate the steady state operating point of a grid-forming VSM QSEM converter without PLL

```

### steadystate_VSM_VC.m

**Signature:** `function conv = steadystate_VSM_VC(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-forming VSM VC converter with PLL

#### Help

```matlab
Calculate the steady state operating point of a grid-forming VSM VC converter with PLL

```

### steadystate_VSM_VC_noPLL.m

**Signature:** `function conv = steadystate_VSM_VC_noPLL(conv,input)`

**Summary:** Calculate the steady state operating point of a grid-forming VSM VC converter without PLL

#### Help

```matlab
Calculate the steady state operating point of a grid-forming VSM VC converter without PLL

```

### steadystate_nonlinear_VSC_GF.m

**Signature:** `function y=steadystate_nonlinear_VSC_GF(x,conv)`

**Summary:** Steady-state nonlinear equations for a grid-following VSC

#### Help

```matlab
Steady-state nonlinear equations for a grid-following VSC

```

### steadystate_nonlinear_VSM_DEM.m

**Signature:** `function y=steadystate_nonlinear_VSM_DEM(x,conv)`

**Summary:** Steady-state nonlinear equations for a grid-forming VSM DEM converter with PLL

#### Help

```matlab
Steady-state nonlinear equations for a grid-forming VSM DEM converter with PLL

```

### steadystate_nonlinear_VSM_DEM_noPLL.m

**Signature:** `function y=steadystate_nonlinear_VSM_DEM_noPLL(x,conv)`

**Summary:** Steady-state nonlinear equations for a grid-forming VSM DEM converter without PLL

#### Help

```matlab
Steady-state nonlinear equations for a grid-forming VSM DEM converter without PLL

```

### steadystate_nonlinear_VSM_QSEM.m

**Signature:** `function y=steadystate_nonlinear_VSM_QSEM(x,conv)`

**Summary:** Steady-state nonlinear equations for a grid-forming VSM QSEM converter with PLL

#### Help

```matlab
Steady-state nonlinear equations for a grid-forming VSM QSEM converter with PLL

```

### steadystate_nonlinear_VSM_QSEM_noPLL.m

**Signature:** `function y=steadystate_nonlinear_VSM_QSEM_noPLL(x,conv)`

**Summary:** Steady-state nonlinear equations for a grid-forming VSM QSEM converter without PLL

#### Help

```matlab
Steady-state nonlinear equations for a grid-forming VSM QSEM converter without PLL

```

### steadystate_nonlinear_VSM_VC.m

**Signature:** `function y=steadystate_nonlinear_VSM_VC(x,conv)`

**Summary:** Steady-state nonlinear equations for a grid-forming VSM VC converter with PLL

#### Help

```matlab
Steady-state nonlinear equations for a grid-forming VSM VC converter with PLL

```

### steadystate_nonlinear_VSM_VC_noPLL.m

**Signature:** `function y=steadystate_nonlinear_VSM_VC_noPLL(x,conv)`

**Summary:** Steady-state nonlinear equations for a grid-forming VSM VC converter without PLL

#### Help

```matlab
Steady-state nonlinear equations for a grid-forming VSM VC converter without PLL

```

