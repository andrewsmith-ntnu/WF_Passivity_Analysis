function [Z_wf,Z_wt,cvtr]=IM_WF(freq,Sb_WF,Sb_WT,n_s,V_wt,P_pu,Q_pu,phi_wt,D_WW,cable_params,tf_params,Z_cable_g)
%% Wind Turbine Impedance

Sb_WT=5e6;   %Nominal power rating of WT
Vdc=1500;
Vnr=V_wt;    %AC voltage, rms
Inr=Sb_WT/sqrt(3)/Vnr;

% P=1;
V1=Vnr*sqrt(2);
Iph=P_pu*Sb_WT/sqrt(3)/Vnr;
I1=Iph*sqrt(2);
phii1=phi_wt;


%% Line and Transformer impedances
length_WW=D_WW;  %km
length_SW=4*D_WW;  %km

r_ww=cable_params(1,1);
l_ww=cable_params(1,2);
c_ww=cable_params(1,3);
r_sw=cable_params(2,1);
l_sw=cable_params(2,2);
c_sw=cable_params(2,3);
R_tf_wt=tf_params(1,1);
L_tf_wt=tf_params(1,2);
R_tf_wf=tf_params(2,1);
L_tf_wf=tf_params(2,2);

%At 33kV
Z_Cf=IM_C(freq,1e-6);     %WT filter capacitance
Z_line_WW=IM_line(freq,r_ww*length_WW,l_ww*length_WW);   %Collection cable series imp
Z_C_WW=IM_C(freq,c_ww*length_WW);      %Collection cable shunt imp
Z_line_SW=IM_line(freq,r_sw*length_SW,l_sw*length_SW);   %Substation-WT cable series imp
Z_C_SW=IM_C(freq,c_sw*length_SW);      %Substation-WT cable shunt imp
Z_trans_WT=IM_line(freq,R_tf_wt,L_tf_wt);    %At 33kV
%At 220kV
Z_trans_WF=IM_line(freq,R_tf_wf,L_tf_wf); %At 220kV
% Z_trans_WF=(2*Z_trans_WF^-1)^-1;    %2 in parallel

%Convert to 33kV
% Z_trans_WF=Z_trans_WF/(33e3/220e3)^2;  


% %% Assemble WF from components, Sainz
% n_t=round(Sb_WF/Sb_WT);
% 
% n_p=round(n_t/n_s);  %Number of strings in parallel
% 
% % %SEQ, grid-following
% % [Z_wt,~]=IM_WT_noCouple(freq,Vdc,V1,I1,phii1);
% 
% %VSM, dq
% [Z_wt,~,cvtr]=IM_WT_VSM(freq,Vph,Iph,P_pu,Q_pu);
% 
% Z_wt=Z_wt*cvtr.pu.Zb*(33e3/V_wt)^2;
% 
% % Z_WT_F=(Z_wt^-1+Z_Cf^-1)^-1+Z_trans_WT;
% Z_WT_F=Z_wt+Z_trans_WT;
% Z_WT_F_L=(((Z_WT_F^-1+Z_C_WW^-1)^-1+Z_line_WW)^-1+Z_C_WW^-1)^-1;
% Z_s=Z_WT_F_L;
% for i_n=1:n_s
%     Z_s=(((Z_s^-1+(Z_WT_F)^-1+Z_C_WW^-1)^-1+Z_line_WW)^-1+Z_C_WW^-1)^-1;
% end
% 
% Z_wf=(n_p*Z_s^-1)^-1;
% Z_wf=(((Z_C_SW^-1+Z_wf^-1)^-1+Z_line_SW)^-1+Z_C_SW^-1)^-1;
% Z_wf=Z_wf+Z_trans_WF;
% 
% 
% 
% 
% %% Impedance as seen by 1 WT
% Z_C_G_T=Z_cable_g+Z_trans_WF;
% Z_GS=(((Z_C_SW^-1+Z_C_G_T^-1)^-1+Z_line_SW)^-1+Z_C_SW^-1)^-1;
% 
% Z_s=Z_WT_F_L;
% for i_n=1:n_s
%     Z_s=(((Z_s^-1+(Z_WT_F)^-1+Z_C_WW^-1)^-1+Z_line_WW)^-1+Z_C_WW^-1)^-1;
% end
% 
% Z_s_m1=Z_WT_F_L;
% for i_n=1:(n_s-1)
%     Z_s_m1=(((Z_s_m1^-1+(Z_WT_F)^-1+Z_C_WW^-1)^-1+Z_line_WW)^-1+Z_C_WW^-1)^-1;
% end
% Z_wf_m1=((n_p-1)*Z_s^-1+Z_s_m1^-1+Z_GS^-1)^-1;
% 
% Z_wf_p1=Z_wf_m1+Z_trans_WT;

%% Assemble WF from components, no collector grid
n_t=round(Sb_WF/Sb_WT);

n_p=round(n_t/n_s);  %Number of strings in parallel

%VSM
[Z_wt,~,cvtr]=IM_WT_VSM(freq,Vnr,Inr,P_pu,Q_pu);

Z_wt=Z_wt*cvtr.pu.Zb*(33e3/V_wt)^2;

Z_WT_F=Z_wt+Z_trans_WT;
Z_s=Z_WT_F;

Z_wf=(n_p*n_s*Z_s^-1)^-1;
Z_wf=Z_wf+Z_trans_WF;

% %% Impedance as seen by 1 WT, no collector
% Z_C_G_T=Z_cable_g+Z_trans_WF;
% Z_GS=(((Z_C_SW^-1+Z_C_G_T^-1)^-1+Z_line_SW)^-1+Z_C_SW^-1)^-1;
% 
% Z_s=Z_WT_F_L;
% for i_n=1:n_s
%     Z_s=(((Z_s^-1+(Z_WT_F)^-1+Z_C_WW^-1)^-1+Z_line_WW)^-1+Z_C_WW^-1)^-1;
% end
% 
% Z_s_m1=Z_WT_F_L;
% for i_n=1:(n_s-1)
%     Z_s_m1=(((Z_s_m1^-1+(Z_WT_F)^-1+Z_C_WW^-1)^-1+Z_line_WW)^-1+Z_C_WW^-1)^-1;
% end
% Z_wf_m1=((n_p-1)*Z_s^-1+Z_s_m1^-1+Z_GS^-1)^-1;
% 
% Z_wf_p1=Z_wf_m1+Z_trans_WT;


%%
% opts_bd=bodeoptions;
% opts_bd.FreqUnits='Hz';
% opts_bd.Grid='on';
% % 
% bodeplot((Z_WT_F^-1*2)^-1,opts_bd);
% % % bodeplot(Z_wf_cable_g/(33e3/220e3)^2,opts_bd);
% hold on;
% % % bodeplot(Z_wt(2,2),opts);
% % % bodeplot(Z_wt/(33e3/V_wt)^2,opts_bd);
% % % hold on;
% % % % bodeplot(Z_line_WW,opts)
% % % bodeplot(Z_conv*cvtr.pu.Zb,opts_bd,'.');
% bodeplot(Z_conv,opts_bd,'- .');