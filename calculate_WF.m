function WF = calculate_WF(WF)

WF.Sb_WF=WF.n_WT*WF.WT.Sb_WT;   %VA, Total WF rated power

WF.WT.cvtr.parameters.V_coll=WF.V_coll;

WF.Trans.n_cables=ceil(WF.Sb_WF/WF.Trans.P_cable_max);   %Number of cables required in parallel
WF.Trans.r_c_eff=WF.Trans.r_c/WF.Trans.n_cables;   %Effective cable resistance/km
WF.Trans.l_c_eff=WF.Trans.l_c/WF.Trans.n_cables;   %Effective cable inductance/km
WF.Trans.c_c_eff=WF.Trans.c_c*WF.Trans.n_cables;   %Effective cable capacitance/km


%% Calculate impedance models
WF.Imp.Z_b_G=(WF.V_trans*sqrt(1/3))^2/WF.Sb_WF;
WF.Imp.Z_b_WT=(WF.V_coll*sqrt(1/3))^2/WF.WT.Sb_WT;
WF.Imp.Z_b_WF=(WF.V_coll*sqrt(1/3))^2/WF.Sb_WF;

%Equivalent grid impedance
WF.Imp.Z_line_G=IM_RL(WF.freq,WF.f1,WF.Trans.rg*WF.Imp.Z_b_G,WF.Trans.lg*(WF.Imp.Z_b_G/(2*pi*WF.f1)));    %Grid equiv.

if strcmp(WF.location,'onshore')
    WF.Imp.Z_line_OH=IM_RL(WF.freq,WF.f1,WF.Trans.r_OH*WF.Imp.Z_b_G,WF.Trans.l_OH*((WF.Imp.Z_b_G/(2*pi*WF.f1))));
    
    %Find equivalent grid impedance as seen by the converter
    Z_eq_on=WF.Imp.Z_line_G+WF.Imp.Z_line_OH;
    Z_f1_on=Z_eq_on.ResponseData(:,:,1);
    r_g_agg_on=(WF.WT.r_tf+WF.r_tf+real(Z_f1_on(1,1)/WF.Imp.Z_b_G)/3);
    l_g_agg_on=(WF.WT.l_tf+WF.l_tf+imag(Z_f1_on(1,1)/WF.Imp.Z_b_G)/3);
    
    %Calculate aggregated WT and WF impedances
    [WF.WT.cvtr_agg,WF.Imp.Z_cvtr_agg]=impedance_WT_agg(WF.freq,WF.WT.cvtr,WF.Sb_WF,r_g_agg_on,l_g_agg_on);
    WF.Imp.Z_tf_WT=IM_RL(WF.Imp.Z_cvtr_agg.Frequency',WF.f1,WF.WT.r_tf*WF.Imp.Z_b_WF,WF.WT.l_tf*(WF.Imp.Z_b_WF/(2*pi*WF.f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    WF.Imp.Z_tf_WF=IM_RL(WF.Imp.Z_cvtr_agg.Frequency',WF.f1,WF.r_tf*WF.Imp.Z_b_WF,WF.l_tf*(WF.Imp.Z_b_WF/(2*pi*WF.f1)))*3;
    WF.Imp.Z_wf=WF.Imp.Z_cvtr_agg+WF.Imp.Z_tf_WT+WF.Imp.Z_tf_WF;

    %Convert to transmission voltage
    WF.Imp.Z_wf=WF.Imp.Z_wf*(WF.V_trans/WF.V_coll)^2;
    %Combine with transmission and grid
    WF.Imp.Z_wf_grid=(WF.Imp.Z_wf^-1+(WF.Imp.Z_line_OH+WF.Imp.Z_line_G)^-1)^-1;
elseif strcmp(WF.location,'offshore')
    WF.Imp.Z_line_C=IM_RL(WF.freq,WF.f1,WF.Trans.r_c_eff*WF.Trans.D_C,WF.Trans.l_c_eff*WF.Trans.D_C);    %Cable series imp
    WF.Imp.Z_C_C=IM_C(WF.freq,WF.f1,WF.Trans.c_c_eff/2*WF.Trans.D_C);                  %Cable shunt imp

    % ####  Reactive compensation, assuming shunt reactor at both ends ####
    WF.l_comp=(1/(WF.Trans.c_c_eff/2*WF.Trans.D_C*2*pi*WF.f1))/WF.Imp.Z_b_G-(WF.l_tf+(WF.Trans.l_c_eff*WF.Trans.D_C*2*pi*WF.f1)/WF.Imp.Z_b_G);
    WF.L_comp=sqrt(2)*WF.l_comp*(WF.Imp.Z_b_G/(2*pi*WF.f1));
    WF.R_comp=0.1*WF.l_comp*WF.Imp.Z_b_G;
    WF.Imp.Z_comp = IM_RL(WF.freq,WF.f1,WF.R_comp,WF.L_comp); 
    r_diss=.1;     %Temporary R to dissipate inrush currents, Simulink

    Z_eq_off=(((WF.Imp.Z_comp^-1+WF.Imp.Z_C_C^-1+WF.Imp.Z_line_G^-1)^-1+WF.Imp.Z_line_C)^-1+WF.Imp.Z_C_C^-1+WF.Imp.Z_comp^-1)^-1;
    Z_f1_on=Z_eq_off.ResponseData(:,:,1);
    r_g_agg_on=(WF.WT.r_tf+WF.r_tf+real(Z_f1_on(1,1)/WF.Imp.Z_b_G)/3);
    l_g_agg_on=(WF.WT.l_tf+WF.l_tf+imag(Z_f1_on(1,1)/WF.Imp.Z_b_G)/3);
    [WF.WT.cvtr_agg,WF.Imp.Z_cvtr_agg]=impedance_WT_agg(WF.freq,WF.WT.cvtr,WF.Sb_WF,r_g_agg_on,l_g_agg_on);
    WF.Imp.Z_tf_WT=IM_RL(WF.Imp.Z_cvtr_agg.Frequency',WF.f1,WF.WT.r_tf*WF.Imp.Z_b_WF,WF.WT.l_tf*(WF.Imp.Z_b_WF/(2*pi*WF.f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
    WF.Imp.Z_tf_WF=IM_RL(WF.Imp.Z_cvtr_agg.Frequency',WF.f1,WF.r_tf*WF.Imp.Z_b_WF,WF.l_tf*(WF.Imp.Z_b_WF/(2*pi*WF.f1)))*3;
    WF.Imp.Z_wf=WF.Imp.Z_cvtr_agg+WF.Imp.Z_tf_WT+WF.Imp.Z_tf_WF;

    %Convert to transmission voltage
    WF.Imp.Z_wf=WF.Imp.Z_wf*(WF.V_trans/WF.V_coll)^2;
    %Combine with transmission and grid
    WF.Imp.Z_wf_grid=(WF.Imp.Z_wf^-1+((WF.Imp.Z_C_C^-1+WF.Imp.Z_line_G^-1+WF.Imp.Z_comp^-1)^-1+WF.Imp.Z_line_C)^-1+WF.Imp.Z_C_C^-1+WF.Imp.Z_comp^-1)^-1;
else
    warning('Location must be onshore or offshore');
end



end
