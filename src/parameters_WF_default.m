% This script initializes default parameters for wind farm (WF) modeling. Use this to change the default values. 
set(0,'defaulttextinterpreter','latex');

% ###### Initialize parameters ######

% ## WT and WF parameters ##
WF.freq=logspace(0,4,800); %frequency range for impedance models 
WF.f1=50;              %Nominal frequency (Hz)
WF.WT.Sb_WT=5e6;       %Wind turbine rated power (VA)
WF.n_WT=200;           %Number of WTs in WF
WF.WT.V_wt=900;        %WT nominal voltage (V)
WF.V_coll=33e3;        %WF collection grid voltage (V)
WF.V_trans=220e3;      %Transmission voltage (V)


%Set wind turbine control parameters
% V_wt=WF.WT.V_wt;
WF.WT.cvtr=initialize_WT_default_gfl(WF);   %Converter parameters, default in GFL mode

% ## Cable and grid impedance ##
WF.Trans.D_C=50;         %Cable transmission distance (km)
%Values for 220kV, 360MW cable (SI)
WF.Trans.r_c=27e-3;      %Cable resistance (Ohm/km)
WF.Trans.l_c=0.386e-3;   %Cable inductance (H/km)
WF.Trans.c_c=0.177e-6;   %Cable capacitance (F/km)

%Overhead line (pu)
WF.Trans.r_OH=0.01;     %Line resistance (pu)
WF.Trans.l_OH=0.1;      %Line inductance (pu)
    
%Grid
WF.Trans.rg=0.01;       %Equivalent grid resistance (pu)
WF.Trans.lg=0.1;        %Equivalent grid inductance (pu)

%Number of transmission cables
WF.Trans.P_cable_max=359e6;  %Max power of single cable (W)

%Transformer values
WF.WT.r_tf=0.0056;      %Wind turbine transformer resistance (pu)
WF.WT.l_tf=0.0392;      %Wind turbine transformer inductance (pu)
WF.r_tf=0.005;          %Wind farm transformer resistance (pu)
WF.l_tf=0.05;           %Wind farm transformer inductance (pu)

%Onshore or offshore
WF.location='onshore';
