function [cvtr_agg] = initialize_WT_agg(cvtr,Sb_WF,r_g_agg,l_g_agg)
% Initializes the aggregated wind turbine converter model based on the single wind turbine converter parameters and aggregated grid-side impedance.

cvtr_agg=cvtr;


%% Electrical ratings
cvtr_agg.parameters.Inr = Sb_WF/cvtr_agg.parameters.Vnr/sqrt(3);

cvtr_agg = base_values(cvtr_agg);

% 
cvtr_agg.parameters.rg=r_g_agg;
cvtr_agg.parameters.lg=l_g_agg;

% Parameters in absolute values
cvtr_agg.parameters.Lf = cvtr_agg.parameters.lf*cvtr_agg.pu.Lb;
cvtr_agg.parameters.Rf = cvtr_agg.parameters.rf*cvtr_agg.pu.Zb;
cvtr_agg.parameters.Cf = cvtr_agg.parameters.cf*cvtr_agg.pu.Cb;
cvtr_agg.parameters.Lg = cvtr_agg.parameters.lg*cvtr_agg.pu.Lb;
cvtr_agg.parameters.Rg = cvtr_agg.parameters.rg*cvtr_agg.pu.Zb;

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
        
        %% Calculate the matrixes
        cvtr_agg = matrix_VSM_DEM_noPLL(cvtr_agg,'VSM');
    case 1
        selector_current_reference = 0;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_DEM(cvtr_agg);
        
        %% Calculate the matrixes
        cvtr_agg = matrix_VSM_DEM(cvtr_agg,'VSM');
    case 2
        selector_current_reference = 1;
        selector_omega_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_QSEM_noPLL(cvtr_agg);
        
        %% Calculate the matrixes
        cvtr_agg = matrix_VSM_QSEM_noPLL(cvtr_agg,'VSM');
    case 3
        selector_current_reference = 1;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_QSEM(cvtr_agg);
        
        %% Calculate the matrixes
        cvtr_agg = matrix_VSM_QSEM(cvtr_agg,'VSM');
    case 4
        selector_current_reference = 2;
        selector_omega_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_VC_noPLL(cvtr_agg);
        
        %% Calculate the matrixes
        cvtr_agg = matrix_VSM_VC_noPLL(cvtr_agg,'VSM');
    case 5
        selector_current_reference = 2;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = steadystate_VSM_VC(cvtr_agg);
        
        %% Calculate the matrixes
        cvtr_agg = matrix_VSM_VC(cvtr_agg,'VSM');
    case 6
        selector_current_reference = 3;
        selector_omega_reference = 1;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr_agg = NR_steadystate_VSC_GF(cvtr_agg);
        
        %% Calculate the matrixes
        cvtr_agg = matrix_VSC_GF(cvtr_agg,'VSM');
end

eig_A = eig(cvtr_agg.A);

warning on
if any(real(eig_A)>0)
    warning('System not stable');
end
