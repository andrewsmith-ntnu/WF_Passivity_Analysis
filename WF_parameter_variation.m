% Main script to run parameter variation cases for WF impedance modeling. 

clear
run('parameters_WF_default.m');
mdl='WF_offshore.slx';

% Calculate impedance models
WF=calculate_WF(WF);


%% Plot nominal passivity of WF and total impedance with grid 
Z1_on=WF.Imp.Z_wf;
Z2_on=WF.Imp.Z_wf_grid;

f=figure;
plot_pass_imp(f,Z1_on,Z2_on,1,1);

%%

WF.location='offshore';
WF=calculate_WF(WF);

Z1_off=WF.Imp.Z_wf;
Z2_off=WF.Imp.Z_wf_grid;

f=figure;
plot_pass_imp(f,Z1_off,Z2_off,1,1);

%% Compare with different grid strengths
run('parameters_WF_default.m');

% #### Define variables ####
WF.WT.cvtr=initialize_WT_default_gfm(WF);
WF.location='offshore';
WF.Trans.D_C=150;

rg1=0.01;
lg1=1;
rg2=0.01;
lg2=2;


% #### Calculate using different grids ####
%First grid strength
WF.Trans.rg=rg1;
WF.Trans.lg=lg1;
WF=calculate_WF(WF);
Z1=WF.Imp.Z_wf;
Z2=WF.Imp.Z_wf_grid;

% cvtr_WF=WF.WT.cvtr_agg;      %Set Simulink converter

f=figure;
plot_pass_imp(f,Z1,Z2,1,1);

%Second grid strength
WF.Trans.rg=rg2;
WF.Trans.lg=lg2;
WF=calculate_WF(WF);
Z1=WF.Imp.Z_wf;
Z2=WF.Imp.Z_wf_grid;
plot_pass_imp(f,Z1,Z2,1,1);
legend('$L_g=1$ pu','','$L_g=2$ pu','interpreter','latex','Location','northeast');


%% Plot GFL vs GFM, onshore
run('parameters_WF_default.m');
f_on=figure;

%GFL
WF.WT.cvtr=initialize_WT_default_gfl(WF);
WF_gfl=calculate_WF(WF);

%GFM
WF.WT.cvtr=initialize_WT_default_gfm(WF);
WF_gfm=calculate_WF(WF);

%Plot
plot_pass_imp(f_on,WF_gfl.Imp.Z_wf,WF_gfl.Imp.Z_wf_grid,1,1);hold on
plot_pass_imp(f_on,WF_gfm.Imp.Z_wf,WF_gfm.Imp.Z_wf_grid,1,1);

legend('GFL','','GFM','interpreter','latex','Location','northeast');

%% Plot GFL vs GFM, offshore
run('parameters_WF_default.m');
f_off=figure;
WF.location='offshore';

%GFL
WF.WT.cvtr=initialize_WT_default_gfl(WF);
WF_gfl=calculate_WF(WF);

%GFM
WF.WT.cvtr=initialize_WT_default_gfm(WF);
WF_gfm=calculate_WF(WF);

%Plot
plot_pass_imp(f_off,WF_gfl.Imp.Z_wf,WF_gfl.Imp.Z_wf_grid,1,1);hold on
plot_pass_imp(f_off,WF_gfm.Imp.Z_wf,WF_gfm.Imp.Z_wf_grid,1,1);

legend('GFL','','GFM','interpreter','latex','Location','northeast');

%% Plot WF and cable impedance with negative zone, vary cable distance
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define distances ####
D_C=10:4:150;

% #### Iterate through each distance ####
p_vary.name='Trans.D_C';
p_vary.val=10:4:150;

WF_on=WF;
WF_on.location='onshore';

WF_off=WF;
WF_off.location='offshore';


% WF=setfield(WF,p_vary.name,60);
for i=1:length(p_vary.val)
    
    %Onshore case
    % WF_on.Trans.D_C=p_vary.val(i);
    WF_on=setNestedField(WF_on,p_vary.name,p_vary.val(i));
    WF_on=calculate_WF(WF_on);
    Z1_array_on{i}=WF_on.Imp.Z_wf;
    Z2_array_on{i}=WF_on.Imp.Z_wf_grid;
    
    %Offshore case
    WF_off=setNestedField(WF_off,p_vary.name,p_vary.val(i));
    WF_off=calculate_WF(WF_off);
    Z1_array_off{i}=WF_off.Imp.Z_wf;
    Z2_array_off{i}=WF_off.Imp.Z_wf_grid;

end

% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];

opts.plot_pass=0;
plot_p_vary(f,{WF_off.Imp.Z_wf},Z2_array_off,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');%,'Ticks',p_vary);
c.Layout.Tile = 'east';
clim([p_vary.val(1) p_vary.val(end)]);
title(c,'Distance (km)','Interpreter','latex');


%% Plot WF and cable impedance with negative zone, vary # turbines
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define number of turbines ####
n_all=(WF.Sb_WF/WF.WT.Sb_WT);
n_t=(n_all/2:4:n_all);

% #### Iterate ####
p_vary.name='n_WT';
p_vary.val=n_t;


%Onshore
WF.location='onshore';
WF=calculate_WF(WF);
WF_partial=WF;
for i=1:length(p_vary.val)
    WF_partial=setNestedField(WF_partial,p_vary.name,p_vary.val(i));
    WF_partial=calculate_WF(WF_partial);
    WF.Imp.Z_wf=WF_partial.Imp.Z_cvtr_agg+WF_partial.Imp.Z_tf_WT+WF.Imp.Z_tf_WF;
    WF.Imp.Z_wf=WF.Imp.Z_wf*(WF.V_trans/WF.V_coll)^2;
    WF.Imp.Z_wf_grid=(WF.Imp.Z_wf^-1+(WF.Imp.Z_line_OH+WF.Imp.Z_line_G)^-1)^-1;
    Z1_array_on{i}=WF.Imp.Z_wf;
    Z2_array_on{i}=WF.Imp.Z_wf_grid;
end

%Onshore
WF.location='offshore';
WF=calculate_WF(WF);
WF_partial=WF;
for i=1:length(p_vary.val)
    WF_partial=setNestedField(WF_partial,p_vary.name,p_vary.val(i));
    WF_partial=calculate_WF(WF_partial);
    WF.Imp.Z_wf=WF_partial.Imp.Z_cvtr_agg+WF_partial.Imp.Z_tf_WT+WF.Imp.Z_tf_WF;
    WF.Imp.Z_wf=WF.Imp.Z_wf*(WF.V_trans/WF.V_coll)^2;
    WF.Imp.Z_wf_grid=(WF.Imp.Z_wf^-1+((WF.Imp.Z_C_C^-1+WF.Imp.Z_line_G^-1+WF.Imp.Z_comp^-1)^-1+WF.Imp.Z_line_C)^-1+WF.Imp.Z_C_C^-1+WF.Imp.Z_comp^-1)^-1;
    Z1_array_off{i}=WF.Imp.Z_wf;
    Z2_array_off{i}=WF.Imp.Z_wf_grid;
end

% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,Z1_array_on,Z2_array_on,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary.val(1)/n_all*100 p_vary.val(end)/n_all*100]);
title(c,'\# Turbines (\%)','interpreter','latex');


%% WF impedance with VSM control parameter variation, Active P ref
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define power references ####
P_pu=0.1:0.1:1;
Q_pu=0;


% #### Iterate ####
p_vary.name='WT.cvtr.ss.pref0';
p_vary.val=P_pu;
WF.WT.cvtr.ss.qref0=Q_pu;

%Onshore case
for i=1:length(p_vary.val)
    WF=setNestedField(WF,p_vary.name,p_vary.val(i));
    WF=calculate_WF(WF);
    Z1_array_on{i}=WF.Imp.Z_wf;
    Z2_array_on{i}=WF.Imp.Z_wf_grid;
end


%Offshore case
WF.location='offshore';
for i=1:length(p_vary.val)
    WF=setNestedField(WF,p_vary.name,p_vary.val(i));
    WF=calculate_WF(WF);
    Z1_array_off{i}=WF.Imp.Z_wf;
    Z2_array_off{i}=WF.Imp.Z_wf_grid;
end


% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,Z1_array_on,Z2_array_on,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary.val(1) p_vary.val(end)]);
title(c,'$P_{ref}$','interpreter','latex');


%% WF impedance with VSM control parameter variation, Reactive Q ref
run('parameters_WF_default.m');
clear Z1_array_on Z1_array_off Z2_array_on Z2_array_off

% #### Define power references ####
P_pu=0.5;
Q_pu=0:0.1:1;


% #### Iterate ####
p_vary.name='WT.cvtr.ss.qref0';
p_vary.val=Q_pu;
WF.WT.cvtr.ss.pref0=P_pu;

%Onshore case
WF.location='onshore';
for i=1:length(p_vary.val)
    WF=setNestedField(WF,p_vary.name,p_vary.val(i));
    WF=calculate_WF(WF);
    Z1_array_on{i}=WF.Imp.Z_wf;
    Z2_array_on{i}=WF.Imp.Z_wf_grid;
end


%Offshore case
WF.location='offshore';
for i=1:length(p_vary.val)
    WF=setNestedField(WF,p_vary.name,p_vary.val(i));
    WF=calculate_WF(WF);
    Z1_array_off{i}=WF.Imp.Z_wf;
    Z2_array_off{i}=WF.Imp.Z_wf_grid;
end


% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,Z1_array_on,Z2_array_on,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary.val(1) p_vary.val(end)]);
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
WF.WT.cvtr.parameters.rg=rg;
lg=linspace(0.01, 1,12);


% #### Iterate ####
p_vary.name='Trans.lg';
p_vary.val=lg;
for i=1:length(p_vary.val)
    if p_vary.val==rg
        WF.WT.cvtr.parameters.rg=p_vary.val(i);
        rg_temp=p_vary.val(i);
        lg_temp=lg;
    elseif p_vary.val==lg
        WF.WT.cvtr.parameters.lg=p_vary.val(i);    
        rg_temp=rg;
        lg_temp=p_vary.val(i);
    end

    %Onshore case
    WF.location='onshore';
    WF.Trans.rg=rg_temp;
    WF.Trans.lg=lg_temp;
    WF=calculate_WF(WF);
    Z1_array_on{i}=WF.Imp.Z_wf;
    Z2_array_on{i}=WF.Imp.Z_wf_grid;

    %Onshore case
    WF.location='offshore';
    WF.Trans.rg=rg_temp;
    WF.Trans.lg=lg_temp;
    WF=calculate_WF(WF);
    Z1_array_off{i}=WF.Imp.Z_wf;
    Z2_array_off{i}=WF.Imp.Z_wf_grid;
end


% #### Plot Results ####
f=figure;
f.Position=f.Position.*[0,0,0.75*2,1.5];
opts.plot_pass=1;
plot_p_vary(f,{WF.Imp.Z_wf},Z2_array_off,1,1,opts);

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary.val(1) p_vary.val(end)]);
title(c,'$L_g (pu)$','interpreter','latex');



%% Blackbox impedance

run('parameters_WF_default.m');
WF.WT.cvtr=initialize_WT_default_gfm(WF);


ref.pref=1;
ref.qref=0;
ref.vref=1;
ref.omegaref=1;

[cvtr_ss,Z_cvtr]=impedance_WT_blackbox(WF.freq,WF.f1,WF.WT.V_wt,WF.V_coll,WF.WT.Sb_WT,WF.n_WT,ref,0.017,0.157);

