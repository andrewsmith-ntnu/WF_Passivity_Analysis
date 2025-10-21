function [cvtr_agg,Z_cvtr_agg] = impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg)
% impedance_WT_agg  Aggregate single-WT model into equivalent aggregated model
%
% [cvtr_agg, Z_cvtr_agg] = impedance_WT_agg(freq, cvtr, Sb_WF, r_g_agg, l_g_agg)
%
% Inputs
%   freq    - frequency vector (Hz) for impedance calculation
%   cvtr    - base converter model struct
%   Sb_WF   - apparent power of the wind farm (per unit or absolute depending on usage)
%   r_g_agg - equivalent grid resistance for aggregated model
%   l_g_agg - equivalent grid inductance for aggregated model
%
% Outputs
%   cvtr_agg    - converter struct updated for aggregated ratings and parameters
%   Z_cvtr_agg  - aggregated converter impedance (sequence-domain FRD)


cvtr_agg = cvtr;


%% Electrical ratings
% Update rated currents for aggregated wind farm
cvtr_agg.parameters.Inr = Sb_WF / cvtr_agg.parameters.Vnr / sqrt(3);

% initialize base values and per-unit scaling
cvtr_agg = base_values(cvtr_agg);

% set aggregated grid parameters
cvtr_agg.parameters.rg = r_g_agg;
cvtr_agg.parameters.lg = l_g_agg;

% Parameters in absolute values (using pu base from base_values)
cvtr_agg.parameters.Lf = cvtr_agg.parameters.lf * cvtr_agg.pu.Lb;
cvtr_agg.parameters.Rf = cvtr_agg.parameters.rf * cvtr_agg.pu.Zb;
cvtr_agg.parameters.Cf = cvtr_agg.parameters.cf * cvtr_agg.pu.Cb;
cvtr_agg.parameters.Lg = cvtr_agg.parameters.lg * cvtr_agg.pu.Lb;
cvtr_agg.parameters.Rg = cvtr_agg.parameters.rg * cvtr_agg.pu.Zb;

%% Controller parameters


% Tuning current loop PU
cvtr_agg.parameters.Tl = cvtr_agg.parameters.lf/(cvtr_agg.parameters.rf*cvtr_agg.pu.Omegab);
cvtr_agg.parameters.Tic = cvtr_agg.parameters.Tl;
cvtr_agg.parameters.kpc = cvtr_agg.parameters.Tl*cvtr_agg.parameters.rf/(2*cvtr_agg.parameters.T_sw);
cvtr_agg.parameters.kic = cvtr_agg.parameters.kpc/cvtr_agg.parameters.Tic;

% Tuning voltage loop PU, VCVSM
cvtr_agg.parameters.av = 2;
cvtr_agg.parameters.Tc = cvtr_agg.parameters.cf/cvtr_agg.pu.Omegab;
cvtr_agg.parameters.Teqc = 2*cvtr_agg.parameters.T_sw;
cvtr_agg.parameters.Tiv = cvtr_agg.parameters.av^2*cvtr_agg.parameters.Teqc;
cvtr_agg.parameters.kpv = cvtr_agg.parameters.Tc/(cvtr_agg.parameters.av*cvtr_agg.parameters.Teqc);
cvtr_agg.parameters.kiv = cvtr_agg.parameters.kpv/cvtr_agg.parameters.Tiv;


% PQ power loop
cvtr_agg.parameters.ap = 2;
cvtr_agg.parameters.Tfp = cvtr_agg.parameters.cf/cvtr_agg.pu.Omegab;
cvtr_agg.parameters.Teqc = 2*cvtr_agg.parameters.T_sw;
cvtr_agg.parameters.Tip = cvtr_agg.parameters.ap^2*cvtr_agg.parameters.Teqc;
cvtr_agg.parameters.kpp = cvtr_agg.parameters.Tfp/(cvtr_agg.parameters.ap*cvtr_agg.parameters.Teqc);
cvtr_agg.parameters.kip = cvtr_agg.parameters.kpp/cvtr_agg.parameters.Tip;
cvtr_agg.parameters.kpq = cvtr_agg.parameters.kpp;
cvtr_agg.parameters.kiq = cvtr_agg.parameters.kip;



%%
switch cvtr_agg.parameters.model_selector
    case 0
        selector_current_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_DEM_noPLL(cvtr_agg);
        
        %% Calculate the matrices
        cvtr_agg = matrix_VSM_DEM_noPLL(cvtr_agg,'VSM');
    case 1
        selector_current_reference = 0;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_DEM(cvtr_agg);
        
        %% Calculate the matrices
        cvtr_agg = matrix_VSM_DEM(cvtr_agg,'VSM');
    case 2
        selector_current_reference = 1;
        selector_omega_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_QSEM_noPLL(cvtr_agg);
        
        %% Calculate the matrices
        cvtr_agg = matrix_VSM_QSEM_noPLL(cvtr_agg,'VSM');
    case 3
        selector_current_reference = 1;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_QSEM(cvtr_agg);
        
        %% Calculate the matrices
        cvtr_agg = matrix_VSM_QSEM(cvtr_agg,'VSM');
    case 4
        selector_current_reference = 2;
        selector_omega_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_VC_noPLL(cvtr_agg);
        
        %% Calculate the matrices
        cvtr_agg = matrix_VSM_VC_noPLL(cvtr_agg,'VSM');
    case 5
        selector_current_reference = 2;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_VC(cvtr_agg);
        
        %% Calculate the matrices
        cvtr_agg = matrix_VSM_VC(cvtr_agg,'VSM');
    case 6
        selector_current_reference = 3;
        selector_omega_reference = 1;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = NR_steadystate_VSC_GF(cvtr_agg);
        
        %% Calculate the matrices
        cvtr_agg = matrix_VSC_GF(cvtr_agg,'VSM');
end

eig_A = eig(cvtr_agg.A);

warning on
if any(real(eig_A)>0)
    warning('System not stable');
end



%% Calculate impedance model
z_cvtr_agg=IM_WT(freq,cvtr_agg);
Z_cvtr_agg=z_cvtr_agg*cvtr_agg.pu.Zb*(cvtr.parameters.V_coll/cvtr_agg.parameters.Vnr)^2;