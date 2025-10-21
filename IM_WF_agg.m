function [Z_wf,Z_wt_agg]=IM_WF_agg(f1,Z_conv_agg,tf_params,Z_b_WF)
%This function calculates the wind farm impedance at the substation, at the
%substation voltage

%f1: fundamental frequency
%Z_conv_agg: aggregated WT converter impedance at collection voltage
%tf_params: transformer parameters in pu at collection voltage
%Z_b_WF: wind farm base impedance at collection voltage. Note different
%base to match Simulink 

%% Line and Transformer impedances

r_tf_wt=tf_params(1,1);
l_tf_wt=tf_params(1,2);
r_tf_wf=tf_params(2,1);
l_tf_wf=tf_params(2,2);


%At collection grid voltage
Z_trans_WT=IM_RL(Z_conv_agg.Frequency',f1,1*r_tf_wt*Z_b_WF,1*l_tf_wt*(Z_b_WF/(2*pi*f1)))*3;   %Note *3 factor to match Simulink transformer model. Due to different base impedance calculation
Z_trans_WF=IM_RL(Z_conv_agg.Frequency',f1,1*r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*f1)))*3;


%% Assemble WF from components, no collector grid

Z_wt_agg=Z_conv_agg+Z_trans_WT;
Z_wf=Z_wt_agg+Z_trans_WF;