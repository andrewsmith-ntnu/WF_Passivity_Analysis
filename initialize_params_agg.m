%%
freq=logspace(0,4,800);

opts_bd=bodeoptions;
opts_bd.FreqUnits='Hz';
opts_bd.Grid='on';
% 
% bodeplot((Z_wf(1,1)),opts);

f1=50;
Sb_WT=5e6;

Sb_WF=200*Sb_WT;
% Sb_WF=12*Sb_WT;
% Sb_WF=36*Sb_WT;
% n_s=6;
n_s=1;
n_GFM=3;
n_GFL=0;
V_coll=33e3;        %Phase-to-phase RMS
V_trans=220e3;      %Phase-to-phase RMS
V_b_coll=V_coll*sqrt(1/3);
V_b_trans=V_trans*sqrt(1/3);
% V_coll=33e3;
% V_trans=V_coll;

D_WW=1;

% Cable and grid impedance
% Z_b_G=(V_trans/sqrt(3/2))^2/Sb_WF;
Z_b_G=(V_b_trans)^2/Sb_WF;
% Z_b_G=(V_trans)^2/(Sb_WT);
D_C=50;    %km
%Values for 220kV, 360MW cable
r_c=27e-3;      %Ohm/km
l_c=0.386e-3;   %H/km
c_c=0.177e-6;   %F/km
% c_c=0.280e-6;   %F/km
% Rg=0.86;t
% Lg=99.03e-3;
% Rg=0.001/Z_b_G;
% Lg=0.01/((Z_b_G/(2*pi*f1)));
% Rg=0.001*Z_b_G;
% Lg=0.01*((Z_b_G/(2*pi*f1)));
%Overhead line
r_OH=0.02;
l_OH=0.1;

Rg=0.01*Z_b_G;
% Lg=0.1*((Z_b_G/(2*pi*f1)));
% Lg=0.1*((Z_b_G/(2*pi*f1)));
% Lg=0.1*((Z_b_G/(2*pi*f1)));
Lg=0.1*((Z_b_G/(2*pi*f1)));


%Number of transmission cables, Pmax=359MW
n_cables=ceil(Sb_WF/359e6);
% n_cables=1;
r_c_eff=r_c/n_cables;
l_c_eff=l_c/n_cables;
c_c_eff=c_c*n_cables;

%At 220kV
Z_line_C=IM_line(freq,r_c_eff*D_C,l_c_eff*D_C);    %Cable series imp
Z_C_C=IM_C(freq,c_c_eff/2*D_C);                  %Cable shunt imp
Z_line_G=IM_line(freq,Rg,Lg);    %Grid equiv.
Z_line_OH=IM_line(freq,r_OH*Z_b_G,l_OH*((Z_b_G/(2*pi*f1))));


% WF impedance

%Impedance params
r_ww=0.098;
l_ww=0.36e-3;
c_ww=0.23e-6;
r_sw=0.041;
l_sw=0.31e-3;
c_sw=0.34e-6;
cable_params=[r_ww,l_ww,c_ww;r_sw,l_sw,c_sw];

% Z_b_WT=(V_coll/sqrt(3/2))^2/Sb_WT;
% Z_b_WF=(V_coll/sqrt(3/2))^2/Sb_WF;

Z_b_WT=V_b_coll^2/Sb_WT;
Z_b_WF=V_b_coll^2/Sb_WF;

% R_tf_wt=1.46;
% L_tf_wt=31e-3;
% R_tf_wf=0.22;
% L_tf_wf=20.63e-3;
% r_tf_wt=R_tf_wt/Z_b_WT;
% l_tf_wt=L_tf_wt/(Z_b_WT/(2*pi*f1));
% r_tf_wf=R_tf_wf/Z_b_WF;
% l_tf_wf=L_tf_wf/(Z_b_WF/(2*pi*f1));

r_tf_wt=0.0056;
l_tf_wt=0.0392;
% r_tf_wt=0.00056;
% l_tf_wt=0.00392;
R_tf_wt=r_tf_wt*Z_b_WT;
L_tf_wt=l_tf_wt*(Z_b_WT/(2*pi*f1));
r_tf_wf=0.005;
l_tf_wf=0.05;
R_tf_wf=r_tf_wf*Z_b_WF;
L_tf_wf=l_tf_wf*(Z_b_WF/(2*pi*f1));

tf_params=[R_tf_wt,L_tf_wt;R_tf_wf,L_tf_wf];

Z_tf_WF=IM_line(freq,R_tf_wf,L_tf_wf); 


% [Z_wf,Z_wt,cvtr]=IM_WF(freq,Sb_WF,Sb_WT,n_s,V_wt,P_pu,Q_pu,phi_wt,D_WW,cable_params,tf_params,Z_cable_g);
[Z_wf,Z_conv,cvtr,Z_wt]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
% [Z_wf,Z_conv,cvtr,Z_wt]=IM_WF_nt(freq,4,Sb_WT,cvtr,D_WW,cable_params,tf_params,V_coll);
% [Z_wf,Z_wt,cvtr_GFL,cvtr_GFM]=IM_WF3(freq,Sb_WF,Sb_WT,n_GFL,n_GFM,cvtr,D_WW,cable_params,tf_params,V_coll);

% for i_GFL=1:n_GFL
%     assignin('base',strcat('cvtr',num2str(i_GFL)),cvtr_GFL);
% end
% 
% for i_GFM=n_GFL:(n_GFL+n_GFM)
%     assignin('base',strcat('cvtr',num2str(i_GFM)),cvtr_GFM);
% end


% run('initialize_adm_WT.m');
% run('initialize_WT.m');
% Z_eq=(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
% Z_eq=Z_line_G+Z_line_C;
Z_eq=Z_line_G;
[~,n_f] = min(abs(Z_eq.Frequency-f1));
Z_f1=Z_eq.ResponseData(:,:,1);


Ib_G=Sb_WF/V_trans/sqrt(3)*sqrt(2);
Z_b_G2=V_trans*sqrt(2/3)/Ib_G;

% r_g_agg=r_tf_wt+r_tf_wf+(r_c_eff*D_C+Rg)/Z_b_G/2;
% l_g_agg=l_tf_wt+l_tf_wf+(l_c_eff*D_C+Lg)/(Z_b_G/(2*pi*f1))/2;
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);

Z_line_g_agg=IM_line(freq,r_g_agg*Z_b_WF,l_g_agg*((Z_b_WF/(2*pi*f1))));


cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);

[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
% [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_trans_WF=IM_line(freq,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*cvtr.pu.fb))); %At 33k

Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wt=Z_wt*(V_trans/V_coll)^2;
cvtr=cvtr_agg;

%Reactive compensation, assuming shunt reactor at both ends
% l_comp=8; %pu
l_comp=(1/(c_c_eff/2*D_C*2*pi*f1))/Z_b_G-(l_tf_wf+(l_c_eff*D_C*2*pi*f1)/Z_b_G);
% l_comp=(1/(c_c_eff*D_C*2*pi*f1))/Z_b_G-((l_c_eff*D_C*2*pi*f1)/Z_b_G);
L_comp=sqrt(2)*l_comp*(Z_b_G/(2*pi*f1));
% L_comp=l_comp*(Z_b_G/(2*pi*f1));
% L_comp=sqrt(2)*((1/(c_c_eff/2*D_C*2*pi*f1))-(l_tf_wf*(Z_b_G)+(l_c_eff*D_C*2*pi*f1)))/(2*pi*f1);
% L_comp=1/(c_c_eff/2*D_C*(2*pi*f1)^2);
% R_comp=0.06*Z_b_G;
R_comp=0.1*l_comp*Z_b_G;
Z_comp = IM_line(freq,R_comp,L_comp); 
r_diss=.1;     %Temporary R to dissipate inrush currents
