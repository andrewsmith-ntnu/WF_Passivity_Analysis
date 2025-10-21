function [Z_wf,Z_conv,cvtr,Z_wt]=IM_WF_nt(freq,n_t,Sb_WT,cvtr,D_WW,cable_params,tf_params,V_coll)
%% Wind Turbine Impedance

cvtr.parameters.Inr=Sb_WT/sqrt(3)/cvtr.parameters.Vnr;



%% Line and Transformer impedances


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


Z_trans_WT=IM_line(freq,R_tf_wt,L_tf_wt);    %At 33k 
Z_trans_WF=IM_line(freq,R_tf_wf,L_tf_wf); %At 33k
% Z_trans_WF=(2*Z_trans_WF^-1)^-1;    %2 in parallel


%% Assemble WF from components, no collector grid

%VSM
[Z_conv,~,cvtr]=IM_WT_all(freq,cvtr);

Z_conv=Z_conv*cvtr.pu.Zb*(V_coll/cvtr.parameters.Vnr)^2/sqrt(2);

Z_wt=Z_conv+Z_trans_WT;
Z_s=Z_wt;

Z_wf=(n_t*Z_s^-1)^-1;
Z_wf=Z_wf+Z_trans_WF;
