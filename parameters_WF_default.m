% clear;
path_for_save='Figures/';
set(0,'defaulttextinterpreter','latex');

% ###### Initialize parameters ######

% ## WT and WF parameters ##
WF.freq=logspace(0,4,800); %frequency range for impedance models 
WF.f1=50;              %Hz, Nominal frequency
WF.WT.Sb_WT=5e6;          %VA, Wind turbine rated power
WF.n_WT=200;           %Number of WTs in WF
WF.WT.V_wt=900;           %V, WT nominal voltage
WF.V_coll=33e3;        %V, WF collection grid voltage
WF.V_trans=220e3;      %V, Transmission voltage


%Set wind turbine control parameters
% V_wt=WF.WT.V_wt;
WF.WT.cvtr=initialize_WT_default_gfl(WF);

% ## Cable and grid impedance ##
WF.Trans.D_C=50;    %km
%Values for 220kV, 360MW cable (SI)
WF.Trans.r_c=27e-3;      %Ohm/km
WF.Trans.l_c=0.386e-3;   %H/km
WF.Trans.c_c=0.177e-6;   %F/km

%Overhead line (pu)
WF.Trans.r_OH=0.01;
WF.Trans.l_OH=0.1;

%Grid
WF.Trans.rg=0.01;
WF.Trans.lg=0.1;

%Number of transmission cables
WF.Trans.P_cable_max=359e6;  %W, max power of single cable

%Transformer values
WF.WT.r_tf=0.0056;
WF.WT.l_tf=0.0392;
WF.r_tf=0.005;
WF.l_tf=0.05;

%Onshore or offshore
WF.location='onshore';
