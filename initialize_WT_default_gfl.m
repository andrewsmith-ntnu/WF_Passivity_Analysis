function cvtr = initialize_WT_default_gfl(WF)
% CC_DEM_noPLL = 0, CC_DEM = 1, CC_QSEM_noPLL = 2, CC_QSEM =3, VCVSM_noPLL, = 4, VCVSM = 5, GF-PQ = 6

cvtr.parameters.model_selector=6;

%% preinitialize variables that may be missing
cvtr.ss.Kappa0 = 0;
cvtr.ss.Xi0 = 0;
cvtr.ss.Epsilonpll0 = 0;
cvtr.ss.DeltaThetapll0 = 0;
cvtr.ss.vpll0 = 0;
cvtr.ss.Xi0 = 0;
cvtr.ss.vomq0 = 0;
cvtr.ss.vomd0 = 0;
cvtr.ss.ilqref0=0;
cvtr.ss.ildref0=0;
cvtr.ss.Omegavsm0=0;
cvtr.ss.GammaSd0=0;
cvtr.ss.GammaSq0=0;

%% Electrical ratings
cvtr.parameters.Vnr = WF.WT.V_wt; 
cvtr.parameters.Vnr_pk = cvtr.parameters.Vnr*sqrt(2)/sqrt(3);
cvtr.parameters.Inr = WF.WT.Sb_WT/cvtr.parameters.Vnr/sqrt(3);
cvtr.parameters.fnr = WF.f1;
cvtr.parameters.Omegan = 1;
cvtr.parameters.Omegag = 1;

cvtr = base_values(cvtr);

% Filter values
cvtr.parameters.lf = .1;
cvtr.parameters.rf = .005;
cvtr.parameters.cf = 3e-2;

%Grid values
cvtr.parameters.rg = 0.01;
cvtr.parameters.lg = 0.1;

% Parameters in absolute values
cvtr.parameters.Lf = cvtr.parameters.lf*cvtr.pu.Lb;
cvtr.parameters.Rf = cvtr.parameters.rf*cvtr.pu.Zb;
cvtr.parameters.Cf = cvtr.parameters.cf*cvtr.pu.Cb;
cvtr.parameters.Lg = cvtr.parameters.lg*cvtr.pu.Lb;
cvtr.parameters.Rg = cvtr.parameters.rg*cvtr.pu.Zb;

%% Controller parameters
cvtr.parameters.f_sw = 10e3; % switching frequency
cvtr.parameters.T_sw = 1/(cvtr.parameters.f_sw);
cvtr.parameters.T_cont = 1/cvtr.parameters.f_sw;     %Control sampling time
cvtr.parameters.Ts_meas = 1/cvtr.parameters.f_sw;
cvtr.parameters.Ts = 1e-5;        %Simulation time step

opts.inp.ts=cvtr.parameters.T_cont;
cvtr.parameters.discrete = 0;     %Discrete or continuous controller

if cvtr.parameters.discrete==0
    cvtr.parameters.rt=cvtr.parameters.Ts;
else
    cvtr.parameters.rt=cvtr.parameters.T_cont;
end


% Tuning current loop PU
cvtr.parameters.Tl = cvtr.parameters.lf/(cvtr.parameters.rf*cvtr.pu.Omegab);
cvtr.parameters.Tic = cvtr.parameters.Tl;
cvtr.parameters.kpc = cvtr.parameters.Tl*cvtr.parameters.rf/(2*cvtr.parameters.T_sw);
% cvtr.parameters.kpc = cvtr.parameters.Tl*cvtr.parameters.rf/(cvtr.parameters.T_sw);
cvtr.parameters.kic = cvtr.parameters.kpc/cvtr.parameters.Tic;

% Tuning voltage loop PU, VCVSM
cvtr.parameters.av = 2;
cvtr.parameters.Tc = cvtr.parameters.cf/cvtr.pu.Omegab;
cvtr.parameters.Teqc = 2*cvtr.parameters.T_sw;
cvtr.parameters.Tiv = cvtr.parameters.av^2*cvtr.parameters.Teqc;
cvtr.parameters.kpv = cvtr.parameters.Tc/(cvtr.parameters.av*cvtr.parameters.Teqc);
cvtr.parameters.kiv = cvtr.parameters.kpv/cvtr.parameters.Tiv;

% Feedforward terms
cvtr.parameters.kffv = 1;

% Tuning pll
cvtr.parameters.apll = 2;
cvtr.parameters.Tfpll = 0.002;
cvtr.parameters.Tipll = cvtr.parameters.apll^2*cvtr.parameters.Tfpll;
cvtr.parameters.kppll = 1/(cvtr.parameters.apll*2*pi*cvtr.parameters.Tfpll);
cvtr.parameters.kppll = cvtr.parameters.kppll/cvtr.pu.Omegab;
cvtr.parameters.kipll = cvtr.parameters.kppll/cvtr.parameters.Tipll;
cvtr.parameters.Omegalppll = 1/cvtr.parameters.Tfpll;

% Active damping
cvtr.parameters.Tfad = 0.005;
cvtr.parameters.Omegaad = 1/cvtr.parameters.Tfad;
cvtr.parameters.kad = 5;

% VSM virtual impedance
cvtr.parameters.rs = .02;
cvtr.parameters.ls = .5;
cvtr.parameters.Omegavo = 500;  %QSEM

% Droop and VSM
cvtr.parameters.ta = 50;
cvtr.parameters.kd = 40;
cvtr.parameters.kdrpOmega = 20;
cvtr.parameters.kdrpq = 0.05;
cvtr.parameters.Omegadf = 50;
cvtr.parameters.Omegaqf = 100;

%Active and reactive power controllers, GFL
cvtr.parameters.ap = 2;
cvtr.parameters.Tfp = cvtr.parameters.cf/cvtr.pu.Omegab;
cvtr.parameters.Teqc = 2*cvtr.parameters.T_sw;
cvtr.parameters.Tip = cvtr.parameters.ap^2*cvtr.parameters.Teqc;
cvtr.parameters.kpp = cvtr.parameters.Tfp/(cvtr.parameters.ap*cvtr.parameters.Teqc);
cvtr.parameters.kip = cvtr.parameters.kpp/cvtr.parameters.Tip;
cvtr.parameters.kpq=cvtr.parameters.kpp;
cvtr.parameters.kiq=cvtr.parameters.kip;


%% References
cvtr.ss.pref0 = 1;                  
cvtr.ss.qref0 = 0;                  
cvtr.ss.vref0 = 1;                  
cvtr.ss.Omegaref0 = 1;              
cvtr.ss.vn0 = cvtr.parameters.Vnr_pk/cvtr.pu.Vb;    
cvtr.ss.vnd0 = cvtr.parameters.Vnr_pk/cvtr.pu.Vb;   
cvtr.ss.vnq0 = 0;                   
cvtr.ss.Omegag0 = 1;                        
cvtr.ss.DeltaThetavsm0 = 0;         

switch cvtr.parameters.model_selector
    case 0
        selector_current_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr = steadystate_VSM_DEM_noPLL(cvtr);
        
        %% Calculate the matrixes
        cvtr = matrix_VSM_DEM_noPLL(cvtr,'VSM');
    case 1
        selector_current_reference = 0;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr = steadystate_VSM_DEM(cvtr);
        
        %% Calculate the matrixes
        cvtr = matrix_VSM_DEM(cvtr,'VSM');
    case 2
        selector_current_reference = 1;
        selector_omega_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr = steadystate_VSM_QSEM_noPLL(cvtr);
        
        %% Calculate the matrixes
        cvtr = matrix_VSM_QSEM_noPLL(cvtr,'VSM');
    case 3
        selector_current_reference = 1;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr = steadystate_VSM_QSEM(cvtr);
        
        %% Calculate the matrixes
        cvtr = matrix_VSM_QSEM(cvtr,'VSM');
    case 4
        selector_current_reference = 2;
        selector_omega_reference = 0;
        enable_pll = 0;
        
        %% Steady state conditions
        cvtr = steadystate_VSM_VC_noPLL(cvtr);
        
        %% Calculate the matrixes
        cvtr = matrix_VSM_VC_noPLL(cvtr,'VSM');
    case 5
        selector_current_reference = 2;
        selector_omega_reference = 0;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr = steadystate_VSM_VC(cvtr);
        
        %% Calculate the matrixes
        cvtr = matrix_VSM_VC(cvtr,'VSM');
    case 6
        selector_current_reference = 3;
        selector_omega_reference = 1;
        enable_pll = 1;
        
        %% Steady state conditions
        cvtr = NR_steadystate_VSC_GF(cvtr);
        
        %% Calculate the matrixes
        cvtr = matrix_VSC_GF(cvtr,'VSM');
end

eig_A = eig(cvtr.A);

if any(real(eig_A)>0)
    warning('System not stable');
end