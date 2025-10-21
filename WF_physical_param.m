% clear;
set(0,'defaulttextinterpreter','latex');

run('parameters_WF_default.m');

%% Plot nominal passivity of WF and total impedance with grid 

Z1=Z_wf;
Z2=Z_wf_cable_grid;
Z2=Z_wf_OH;

f=figure;
plot_pass_imp(f,Z1,Z2,1,1);
%% Compare with different grid strengths
run('parameters_WF_default.m');

mdl='WF_offshore.slx';

% #### Define variables ####
run('initialize_WT_default_gfm.m'); %Use GFM
rg1=0.01;
lg1=1;
rg2=0.01;
lg2=2;


% #### Calculate using different grids ####
Z1=Z_wf;

%First grid strength
Z_line_G=IM_RL(freq,f1,rg1*Z_b_G,lg1*(Z_b_G/(2*pi*f1)));
Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
[cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z2=Z_wf_cable_grid;
cvtr_WF=cvtr_agg;      %Set Simulink converter

f=figure;
plot_pass_imp(f,Z1,Z2,1,1);


%Second grid strength
Z_line_G=IM_RL(freq,f1,rg2*Z_b_G,lg2*(Z_b_G/(2*pi*f1)));
% Z_eq=Z_line_G+Z_line_OH;
Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
[cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z2=Z_wf_cable_grid;

plot_pass_imp(f,Z1,Z2,1,1);
legend('$L_g=1$ pu','','$L_g=2$ pu','interpreter','latex','Location','northeast');



%% Plot GFL vs GFM, onshore
run('parameters_WF_default.m');
f_on=figure;

%Onshore case
run('initialize_WT_default_gfl.m');
cvtr.parameters.model_selector=6;
Z_eq=Z_line_G+Z_line_OH;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
[cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
Z_wf_gfl=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
Z_wf_gfl=Z_wf_gfl*(V_trans/V_coll)^2;
Z_wf=Z_wf_gfl;
Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z_on_gfl=Z_wf_OH;

run('initialize_WT_default_gfm.m');
cvtr.parameters.model_selector=3;
[cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
Z_wf_gfm=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
Z_wf_gfm=Z_wf_gfm*(V_trans/V_coll)^2;
Z_wf=Z_wf_gfm;
Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z_on_gfm=Z_wf_OH;

plot_pass_imp(f_on,Z_wf_gfl,Z_on_gfl,1,1);hold on
plot_pass_imp(f_on,Z_wf_gfm,Z_on_gfm,1,1);

legend('GFL','','GFM','interpreter','latex','Location','northeast');

%% Plot GFL vs GFM, offshore
run('parameters_WF_default.m');
f_off=figure;

run('initialize_WT_default_gfl.m');
cvtr.parameters.model_selector=6;
Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
[cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
Z_wf_gfl=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
Z_wf_gfl=Z_wf_gfl*(V_trans/V_coll)^2;
Z_wf=Z_wf_gfl;
Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_off_gfl=Z_wf_cable_grid;

run('initialize_WT_default_gfm.m');
cvtr.parameters.model_selector=3;
cvtr_agg=initialize_WT_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
z_cvtr_agg=IM_WT(freq,cvtr_agg);
Z_cvtr_agg=z_cvtr_agg*cvtr_agg.pu.Zb*(V_coll/cvtr_agg.parameters.Vnr)^2;
[Z_wf_gfm,~]=IM_WF_agg(f1,Z_cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],Z_b_WF);
Z_wf_gfm=Z_wf_gfm*(V_trans/V_coll)^2;
Z_wf=Z_wf_gfm;
Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_off_gfm=Z_wf_cable_grid;

plot_pass_imp(f_off,Z_wf_gfl,Z_off_gfl,1,1);hold on
plot_pass_imp(f_off,Z_wf_gfm,Z_off_gfm,1,1);

legend('GFL','','GFM','interpreter','latex','Location','northeast');


%% Plot WF and cable impedance with negative zone, vary cable distance
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define distances ####
D_C=10:4:150;

% #### Iterate through each distance ####
p_vary=D_C;
for i=1:length(p_vary)
    Z_line_C=IM_RL(freq,f1,r_c_eff*D_C(i),l_c_eff*D_C(i));    %Cable series imp
    Z_C_C=IM_C(freq,f1,c_c_eff/2*D_C(i));                  %Cable shunt imp
    l_comp=(1/(c_c_eff/2*D_C(i)*2*pi*f1))/Z_b_G-(l_tf_wf+(l_c_eff*D_C(i)*2*pi*f1)/Z_b_G);
    L_comp=sqrt(2)*l_comp*(Z_b_G/(2*pi*f1));
    R_comp=0.1*l_comp*Z_b_G;
    Z_comp = IM_RL(freq,f1,R_comp,L_comp); 
    Z_line_OH=IM_RL(freq,f1,r_OH*Z_b_G*D_C(i)/50,l_OH*((Z_b_G/(2*pi*f1)))*D_C(i)/50);

    %Onshore case
    Z_eq=Z_line_G+Z_line_OH;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
    Z1_array_on{i}=Z_wf;
    Z2_array_on{i}=Z_wf_OH;
    
    %Offshore case
    Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z1_array_off{i}=Z_wf;
    Z2_array_off{i}=Z_wf_cable_grid;

end

% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];

opts.plot_pass=0;
plot_p_vary(f,{Z_wf},Z2_array_off,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');%,'Ticks',p_vary);
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'Distance (km)','Interpreter','latex');


%% Plot WF and cable impedance with negative zone, vary # turbines
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off
Z=Z_wf;

% #### Define number of turbines ####
n_all=(Sb_WF/Sb_WT);
n_t=(n_all/2:4:n_all);


% #### Iterate ####
p_vary=n_t;

%Onshore
Z_eq=Z_line_G+Z_line_OH;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg_on=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg_on=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);

%Offshore
Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg_off=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg_off=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
for i=1:length(p_vary)
    Sb_partial=p_vary(i)/n_all*Sb_WF;

    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_partial,r_g_agg_on,l_g_agg_on);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF*Sb_WF/Sb_partial,l_tf_wt*(Sb_WF/Sb_partial/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
    Z1_array_on{i}=Z_wf;
    Z2_array_on{i}=Z_wf_OH;

    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_partial,r_g_agg_off,l_g_agg_off);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF*Sb_WF/Sb_partial,l_tf_wt*(Z_b_WF*Sb_WF/Sb_partial/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
    Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z1_array_off{i}=Z_wf;
    Z2_array_off{i}=Z_wf_cable_grid;
   
end

% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,Z1_array_on,Z2_array_on,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1)/n_all*100 p_vary(end)/n_all*100]);
title(c,'\# Turbines (\%)','interpreter','latex');


%% WF impedance with VSM control parameter variation, Active P ref
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define power references ####
P_pu=0.1:0.1:1;
Q_pu=0;


% #### Iterate ####
p_vary=P_pu;
cvtr.ss.qref0=Q_pu;

%Onshore case
Z_eq=Z_line_G+Z_line_OH;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);

for i=1:length(p_vary)
    cvtr.ss.pref0=p_vary(i);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
    Z1_array_on{i}=Z_wf;
    Z2_array_on{i}=Z_wf_OH;
end


%Offshore case
Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);

for i=1:length(p_vary)
    cvtr.ss.pref0=p_vary(i);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z1_array_off{i}=Z_wf;
    Z2_array_off{i}=Z_wf_cable_grid;
end


% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,Z1_array_on,Z2_array_on,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$P_{ref}$','interpreter','latex');


%% WF impedance with VSM control parameter variation, Reactive Q ref
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define power references ####
P_pu=0.5;
Q_pu=0:0.1:1;


% #### Iterate ####
p_vary=Q_pu;
cvtr.ss.pref0=P_pu;

%Onshore case
Z_eq=Z_line_G+Z_line_OH;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);

for i=1:length(p_vary)
    cvtr.ss.qref0=p_vary(i);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
    Z1_array_on{i}=Z_wf;
    Z2_array_on{i}=Z_wf_OH;
end


%Offshore case
Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);

for i=1:length(p_vary)
    cvtr.ss.qref0=p_vary(i);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z1_array_off{i}=Z_wf;
    Z2_array_off{i}=Z_wf_cable_grid;
end

% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,Z1_array_on,Z2_array_on,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$Q_{ref}$','interpreter','latex');


%% WF impedance with VSM control parameter variation, grid impedance
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define grid range ####
% %Rg
% rg=linspace(0.001,0.1,9);
% lg=0.1;
% cvtr.parameters.lg=lg;

% Lg
rg=0.01;
cvtr.parameters.rg=rg;
lg=linspace(0.01, 1,12);


% #### Iterate ####
p_vary=lg;
for i=1:length(p_vary)
    if p_vary==rg
        cvtr.parameters.rg=p_vary(i);
        rg_temp=p_vary(i);
        lg_temp=lg;
    elseif p_vary==lg
        cvtr.parameters.lg=p_vary(i);    
        rg_temp=rg;
        lg_temp=p_vary(i);
    end
    Z_line_G=IM_RL(freq,f1,rg_temp*Z_b_G,lg_temp*(Z_b_G/(2*pi*f1)));

    %Onshore case
    Z_eq=Z_line_G+Z_line_OH;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
    Z1_array_on{i}=Z_wf;
    Z2_array_on{i}=Z_wf_OH;

    %Offshore case
    Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    [cvtr_agg,Z_cvtr_agg]=impedance_WT_agg(freq,cvtr,Sb_WF,r_g_agg,l_g_agg);
    Z_trans_WT=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wt*Z_b_WF,l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    Z_trans_WF=IM_RL(Z_cvtr_agg.Frequency',f1,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;
    Z_wf=Z_cvtr_agg+Z_trans_WT+Z_trans_WF;
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_cable_grid=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z1_array_off{i}=Z_wf;
    Z2_array_off{i}=Z_wf_cable_grid;
end


% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,{Z_wf},Z2_array_off,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$L_g (pu)$','interpreter','latex');

