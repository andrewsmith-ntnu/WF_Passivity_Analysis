% clear;
path_for_save='Figures/';
set(0,'defaulttextinterpreter','latex');
% run('initialize_WT_default_GFM.m');
% run('initialize_WT_default_GFL.m');
% run('initialize_WT_default_040325.m');
% run('initialize_WT_default_120325.m');
% run('initialize_WT_default_060525.m');
% run('initialize_WT_default_180625.m');
run('initialize_WT_default_gfl.m');

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
V_coll=33e3;
V_trans=220e3;
V_b_coll=33e3*sqrt(1/3);
V_b_trans=220e3*sqrt(1/3);
% V_coll=33e3;
% V_trans=V_coll;

D_WW=1;

% Cable and grid impedance
% Z_b_G=(V_trans/sqrt(3/2))^2/Sb_WF;
Z_b_G=(V_b_trans)^2/Sb_WF;
% Z_b_G=(V_b_trans)^2/(Sb_WF/3);
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
r_OH=0.01;
l_OH=0.1;

Rg=0.01*Z_b_G;
% Rg=0.1*Z_b_G;
Lg=0.1*((Z_b_G/(2*pi*f1)));
% Lg=0.01*((Z_b_G/(2*pi*f1)));
% Lg=0.11*((Z_b_G/(2*pi*f1)));
% Lg=1*((Z_b_G/(2*pi*f1)));
Lg2=Lg;
Lg1=Lg;


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

% Z_b_WT=V_b_coll^2/Sb_WT;
% Z_b_WF=V_b_coll^2/(Sb_WF/3);

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
% [Z_wf_og,Z_conv,cvtr,Z_wt]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
% [Z_wf,Z_conv,cvtr,Z_wt]=IM_WF_nt(freq,4,Sb_WT,cvtr,D_WW,cable_params,tf_params,V_coll);
% [Z_wf,Z_wt,cvtr_GFL,cvtr_GFM]=IM_WF3(freq,Sb_WF,Sb_WT,n_GFL,n_GFM,cvtr,D_WW,cable_params,tf_params,V_coll);
% Z_wf_og=Z_wf_og*(V_trans/V_coll)^2;

% for i_GFL=1:n_GFL
%     assignin('base',strcat('cvtr',num2str(i_GFL)),cvtr_GFL);
% end
% 
% for i_GFM=n_GFL:(n_GFL+n_GFM)
%     assignin('base',strcat('cvtr',num2str(i_GFM)),cvtr_GFM);
% end



%Offshore case
% Z_eq=Z_line_G+Z_line_C;
%Onshore case
Z_eq=Z_line_G+Z_line_OH;
% Z_eq=Z_line_G;

[~,n_f] = min(abs(Z_eq.Frequency-f1));
Z_f1=Z_eq.ResponseData(:,:,1);


% r_g_agg=r_tf_wt+r_tf_wf+(r_c_eff*D_C+Rg)/Z_b_G/2;
% l_g_agg=l_tf_wt+l_tf_wf+(l_c_eff*D_C+Lg)/(Z_b_G/(2*pi*f1))/2;
% r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/(2*pi));
% l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/(2*pi));

% r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
% % l_g_agg=(l_tf_wt+l_tf_wf);%+0.4);
% l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/(3));

r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
% r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G));
% l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G));


cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);

[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_trans_WF=IM_line(freq,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*cvtr.pu.fb))); %At 33k
cvtr1=cvtr_agg;

Z_conv=Z_conv*(V_trans/V_coll)^2;
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wt=Z_wt*(V_trans/V_coll)^2;

Z_shifted=Z_wf;
phase_shift=-2*cvtr_agg.ss.DeltaThetapll0; %rad

for f=1:length(Z_shifted.Frequency)
    % Z_shifted.ResponseData(1,1,f)=Z_shifted.ResponseData(1,2,f)*exp(1i*phase_shift);
    Z_shifted.ResponseData(1,2,f)=Z_shifted.ResponseData(1,2,f)*exp(1i*phase_shift);
    Z_shifted.ResponseData(2,1,f)=Z_shifted.ResponseData(2,1,f)*exp(-1i*phase_shift);
    % Z_shifted.ResponseData(2,2,f)=Z_shifted.ResponseData(2,1,f)*exp(-1i*phase_shift);
end
Z_wf=Z_shifted;

%Reactive compensation, assuming shunt reactor at both ends
% l_comp=8; %pu
l_comp=(1/(c_c_eff/2*D_C*2*pi*f1))/Z_b_G-(l_tf_wf+(l_c_eff*D_C*2*pi*f1)/Z_b_G);
% l_comp=(1/(c_c_eff*D_C*2*pi*f1))/Z_b_G-((l_c_eff*D_C*2*pi*f1)/Z_b_G);
L_comp=sqrt(2)*l_comp*(Z_b_G/(2*pi*f1));
% L_comp=sqrt(2)*((1/(c_c_eff/2*D_C*2*pi*f1))-(l_tf_wf*(Z_b_G)+(l_c_eff*D_C*2*pi*f1)))/(2*pi*f1);
% L_comp=1/(c_c_eff/2*D_C*(2*pi*f1)^2);
% R_comp=0.06*Z_b_G;
R_comp=0.1*l_comp*Z_b_G;
Z_comp = IM_line(freq,R_comp,L_comp); 
r_diss=.1;     %Temporary R to dissipate inrush currents


Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_cable_g=Z_wf_cable+Z_line_G;
Z_wf_cable_g2=Z_wf+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_cable2=Z_wf+((Z_C_C^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
Z_wf_cable_g3_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_cable_g4=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_grid_wf_ser=Z_line_G+(((Z_C_C^-1+Z_wf^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
% Z_wf_cable_g5=(Z_wt^-1+(Z_trans_WF+Z_line_C+Z_line_G)^-1)^-1;
% Z_wf_cable_g5_ser=Z_wf+Z_line_C+Z_line_G;
% Z_cable_wc=(((Z_C_C^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1+Z_tf_WF;
% Z_cable=((Z_C_C+Z_line_C)^-1+Z_C_C^-1)^-1;
% Z_grid_eq=(((Z_wf^-1+Z_C_C^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1+Z_line_G;
% Z_cable_g=((((Z_C_C^-1+Z_line_G^-1)^-1)+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_g=Z_wf+Z_line_C+Z_line_G;
% Z_wt_g=Z_wt+Z_trans_WF+Z_C_C;
% Z_pcc=((((Z_wf^-1+Z_C_C^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1+Z_line_G^-1)^-1;
% Z_wf_cable_g=(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1+Z_wf;

% Z_wf_cable_parall=(Z_wf^-1+Z_cable_g^-1)^-1;
Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z_wf_OH_ser=Z_wf+(Z_line_OH+Z_line_G);
% % % Impedance from a WT
% Z_wf_wt_p=(Z_wf_p1^-1+Z_wt^-1)^-1;
% Z_wf_wt_s=Z_wf_p1+Z_wt;


Z_b_WF_tf=(V_coll/sqrt(3))^2/(Sb_WF/3);
% Z_b_WF=(V_coll/sqrt(3))^2/(Sb_WF);

%At 33kV
Z_trans_WT=IM_line(freq,r_tf_wt*Z_b_WF_tf,l_tf_wt*(Z_b_WF_tf/(2*pi*cvtr.pu.fb)))*(V_trans/V_coll)^2;    %At 33k 
Z_trans_WF=IM_line(freq,r_tf_wf*Z_b_WF_tf,l_tf_wf*(Z_b_WF_tf/(2*pi*cvtr.pu.fb)))*(V_trans/V_coll)^2; %At 33k

Z_conv_wf=(Z_conv^-1+(Z_trans_WT+Z_trans_WF+Z_line_OH+Z_line_G)^-1)^-1;
Z_conv_wf_ser=Z_conv+(Z_trans_WT+Z_trans_WF+Z_line_OH+Z_line_G);

% Z_filt_L=IM_line(freq,cvtr.parameters.Rf,cvtr.parameters.Lf);    %Cable series imp
% Z_filt_C=IM_C(freq,cvtr.parameters.Cf);   
% Z_test=(Z_filt_L^-1+Z_filt_C^-1)^-1*cvtr.pu.Zb*(V_coll/cvtr.parameters.Vnr)^2;
% Z_trans_WT=IM_line(freq,R_tf_wt,L_tf_wt);
% Z_test2=Z_test+Z_trans_WT;

Z=Z_wf;
pass=mimo_passivity(Z);
neg_pass_freq_i=find(pass<0);
neg_pass_freq=Z.Frequency(neg_pass_freq_i(end));

% %% Plot impedance
opts_bd=bodeoptions;
opts_bd.FreqUnits='Hz';
opts_bd.Grid='on';
opts_bd.MagUnits='abs';
opts_bd.MagScale='log';
opts_bd.IOGrouping='none';
% figure;
% Z_wt=Y_wt^-1;
% % 
% figure;
% bodeplot(Z_wf,opts_bd);
% bodeplot(Z_wf_cable_g4_ser,opts_bd);
% % bodeplot(Z_wf_cable_g(1,1),opts_bd);
% % bodeplot(Z_wf_cable_parall(1,1),opts_bd);
% bodeplot(Z2,opts_bd);hold on

% bodeplot(Z2,opts_bd,'*');
% bodeplot(Z_wf_cable_g3_ser,opts_bd);
% hold on
% % bodeplot(Z_wf_cable_g/(33e3/220e3)^2,opts_bd);
% hold on;
% % bodeplot(Z_wf+Z_line_G,opts_bd);
% bodeplot((Z_wt*(V_wt/V_coll)^2),opts_bd);
% hold on;

% bodeplot((Z_wt(1,1)^-1*10)^-1,opts_bd);
% bodeplot((real(Z_C_C)),opts_bd);
% bodeplot(Z_line_C,opts_bd)

% bodeplot(Z_cable_g,opts_bd)
% bodeplot(Z_conv*(V_coll/V_trans)^2,opts_bd,'- .');
%%
figure;
Z2=Z_wf_OH;
% bodeplot(Z1,opts_bd);hold on
bodeplot(Z2,opts_bd);hold on

%
% 
% p=loglog(fout,squeeze(mag1(o,i,:)));
% hold on
% grid on
% p2=loglog(fout,squeeze(mag2(o,i,:)));
% xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)|$','Interpreter','latex');

%% Plot pass and imag part of Z
% Z=Z_wt;
% Z=Z_wf;
Z=Z_wf;
% Z=Z_wf_OH;
% Z2=Z;
% Z2=Z_wf_g;
% Z2=Z_wf_cable_g3_ser;
% Z2=Z_wf_cable_g3;
% Z2=Z_wf_line_g;
% Z2=Z_wf;
% Z2=Z_cable_wc;
% Z2=Z_wf_cable_g4;
% Z2=Z_grid_wf_ser;
% Z=Z_wf_OH_ser;
Z2=Z_wf_OH;
% Z2=Z_conv_wf_ser;


f=figure;
t=tiledlayout(2,1);
% f=gcf;
f.Position=f.Position.*[0,0,1.5,1.5];

% pass=real(Z.ResponseData);
im_part=imag(Z2.ResponseData);

clear mag;
clear pass;
    %MIMO passivity

mag=abs(Z2.ResponseData);
% pass=real(Z(out,in).ResponseData);
pass=mimo_passivity(Z);

o=1;
i=1;
fout=Z.Frequency;
nexttile(t);
p=semilogx(fout,pass);
ylabel('Passivity');
grid on;
p.LineWidth=1.5;
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;

nexttile(t);
% p=semilogx(fout,squeeze(im_part(o,i,:)));
grid on
p=loglog(fout,squeeze(mag(o,i,:)));
grid on
% p3=loglog(fout,squeeze(mag3(out,in,:)));
% xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)|$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$\Im[Z_{pp}(j\omega)]$','Interpreter','latex');
p.LineWidth=1.5;

ax=gca;
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;


line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');

%%
Z2=Z_wf_cable_g4;

f=figure;
% f=gcf;
f.Position=f.Position.*[0,0,1.5,1.5*0.5];

% pass=real(Z.ResponseData);
im_part=imag(Z2.ResponseData);
mag=abs(Z2.ResponseData);


o=1;
i=1;
% p=semilogx(fout,squeeze(im_part(o,i,:)));
grid on
p=loglog(fout,squeeze(mag(o,i,:)));
grid on
% p3=loglog(fout,squeeze(mag3(out,in,:)));
% xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)|$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$\Im[Z_{pp}(j\omega)]$','Interpreter','latex');
p.LineWidth=1.5;

ax=gca;
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;


line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');

% exportgraphics(f,[path_for_save,'Z_offshore_extreme','.pdf']);
%% Compare with different grid
% run('initialize_WT_default_gfl.m');
Z=Z_wf;

Rg=0.01*Z_b_G;
Lg=0.1*((Z_b_G/(2*pi*f1)));
Z_line_G=IM_line(freq,Rg,Lg);   
Z_eq=Z_line_G+Z_line_OH;
% Z_eq=Z_line_G;
% [~,n_f] = min(abs(Z_eq.Frequency-f1));
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3;
l_g_agg=l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3;
cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_OH_ser=Z_wf+(Z_line_OH+Z_line_G);
Z_wf_OH=((Z_wf)^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z2=Z_wf_OH;
Lg1=Lg;
cvtr1=cvtr_agg;

f=figure;
t=tiledlayout(2,1);
t.TileSpacing='compact';
% f=gcf;
f.Position=f.Position.*[0,0,1.5,1.5];

% pass=real(Z.ResponseData);
im_part=imag(Z2.ResponseData);
mag=abs(Z2.ResponseData);
pass=mimo_passivity(Z);
neg_pass_freq_i=find(pass<0);
neg_pass_freq=Z.Frequency(neg_pass_freq_i(end));

o=1;
i=1;
fout=Z.Frequency;
nexttile(t);
p=semilogx(fout,pass);
% p=semilogx(fout,pass,'Color',[0.8500 0.3250 0.0980]);
ylabel('Passivity');
grid on;
p.LineWidth=1.5;
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
set(gca,'XTickLabel',[]);
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;
line([neg_pass_freq neg_pass_freq],[ax.YLim(1) ax.YLim(end)],'Color','k','LineStyle','--');

nexttile(t);
% p=semilogx(fout,squeeze(im_part(o,i,:)),'Color','#22976C');
p=loglog(fout,squeeze(mag(o,i,:)),'Color','#22976C');
grid on
hold on
xlabel('Frequency (Hz)','interpreter','latex');
% ylabel('$\Im[Z_{pp}(j\omega)]$','Interpreter','latex');
ylabel('$|Z_{pp}(j\omega)|$','Interpreter','latex');
p.LineWidth=1.5;

ax=gca;
f_size=16;
ax.FontSize=f_size;



Rg=0.01*Z_b_G;
Lg=0.8*((Z_b_G/(2*pi*f1)));
Z_line_G=IM_line(freq,Rg,Lg);   
Z_eq=Z_line_G+Z_line_OH;
% Z_eq=Z_line_G;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3;
l_g_agg=l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3;
cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_OH_ser=Z_wf+(Z_line_OH+Z_line_G);
Z_wf_OH=((Z_wf)^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z2=Z_wf_OH;
Lg2=Lg;


im_part=imag(Z2.ResponseData);
mag=abs(Z2.ResponseData);
% p2=semilogx(fout,squeeze(im_part(o,i,:)),'Color',[0.9290 0.6940 0.1250]);
p2=loglog(fout,squeeze(mag(o,i,:)),'Color',[0.9290 0.6940 0.1250]);
p2.LineWidth=1.5;
ax=gca;
% ylim([-50 900])
line([neg_pass_freq neg_pass_freq],[ax.YLim(1) ax.YLim(end)],'Color','k','LineStyle','--');
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
legend('$L_g=0.1$ pu','$L_g=0.8$ pu','Passivity boundary','interpreter','latex','Location','southeast');

% exportgraphics(f,[path_for_save,'Z_onshore_case_gfl','.pdf']);
% exportgraphics(f,[path_for_save,'Z_onshore_case_gfm','.pdf']);

%% Plot GFL vs GFM

f=figure;
t=tiledlayout(3,1);
t.TileSpacing='compact';
f.Position=f.Position.*[0,0,1.5,1.3];

%Onshore case
Z_eq=Z_line_G+Z_line_OH;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
run('initialize_WT_default_gfl.m');
cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_OH_ser=Z_wf+Z_line_OH+Z_line_G;
Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z_on_gfl=Z_wf_OH;

run('initialize_WT_default_gfm.m');
cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_OH_ser=Z_wf+Z_line_OH+Z_line_G;
Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z_on_gfm=Z_wf_OH_ser;

%Offshore case
Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_f1=Z_eq.ResponseData(:,:,1);
r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
run('initialize_WT_default_gfl.m');
cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_wf_cable_g4=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_off_gfl=Z_wf_cable_g4;
Z_wf_gfl=Z_wf;

run('initialize_WT_default_gfm.m');
cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_wf_cable_g4=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
Z_off_gfm=Z_wf_cable_g4;
Z_wf_gfm=Z_wf;


% pass=real(Z.ResponseData);
im_part_on_gfl=imag(Z_on_gfl.ResponseData);
im_part_on_gfm=imag(Z_on_gfm.ResponseData);
im_part_off_gfl=imag(Z_off_gfl.ResponseData);
im_part_off_gfm=imag(Z_off_gfm.ResponseData);

mag_on_gfl=abs(Z_on_gfl.ResponseData);
mag_on_gfm=abs(Z_on_gfm.ResponseData);
mag_off_gfl=abs(Z_off_gfl.ResponseData);
mag_off_gfm=abs(Z_off_gfm.ResponseData);


%MIMO passivity
pass_gfl=mimo_passivity(Z_wf_gfl);
pass_gfm=mimo_passivity(Z_wf_gfm);

Z=Z_wf_gfl;
o=1;
i=1;
fout=Z.Frequency;
nexttile(t,1);
p=semilogx(fout,pass_gfl);
ylabel('Passivity','Interpreter','latex');
grid on; hold on
p.LineWidth=1.5;
p2=semilogx(fout,pass_gfm);
p2.LineWidth=1.5;
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
set(gca,'XTickLabel',[]);
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;



nexttile(t,2);
% p=semilogx(fout,squeeze(im_part_on_gfl(o,i,:)));
p=loglog(fout,squeeze(mag_on_gfl(o,i,:)));
p.LineWidth=1.5;
grid on; hold on
% p2=semilogx(fout,squeeze(im_part_on_gfm(o,i,:)));
p2=loglog(fout,squeeze(mag_on_gfm(o,i,:)));
p2.LineWidth=1.5;
% ylabel('$\Im\{Z_{on,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{on,pp}(j\omega)|$','Interpreter','latex');
ax=gca;
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');

nexttile(t,3);
% p=semilogx(fout,squeeze(im_part_off_gfl(o,i,:)));
p=loglog(fout,squeeze(mag_off_gfl(o,i,:)));
p.LineWidth=1.5;
grid on; hold on
% p2=semilogx(fout,squeeze(im_part_off_gfm(o,i,:)));
p2=loglog(fout,squeeze(mag_off_gfm(o,i,:)));
p2.LineWidth=1.5;
xlabel('Frequency (Hz)','interpreter','latex');
% ylabel('$\Im\{Z_{on,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{off,pp}(j\omega)|$','Interpreter','latex');
ax=gca;
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
legend('GFL','GFM','Location','southeast','interpreter','latex');

% exportgraphics(f,[path_for_save,'Z_default','.pdf']);

%% Analytical vs Sim Impedance
%GFL
run('initialize_WT_default_gfl.m');
run('initialize_params_agg.m');
% cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);

load('vals_pn_adm_WF_060525_gfl');

V_pn=vals{1};
I1_pn=vals{3};
sig=vals{5};

clear Y_c Z_c
for f=1:length(sig(1).freqs)
    Z_c(:,:,f)=[V_pn{1,f}, V_pn{2,f}]*([I1_pn{1,f}, I1_pn{2,f}]^-1);
end

Z_conv_exp=frd(-Z_c,sig(1).freqs);
Z_conv_exp.FrequencyUnit='Hz';

Z_shifted=Z_conv_exp;
phase_shift=4; %rad

for f=1:length(Z_shifted.Frequency)
    Z_shifted.ResponseData(1,2,f)=Z_shifted.ResponseData(1,2,f)*exp(1i*phase_shift);
    Z_shifted.ResponseData(2,1,f)=Z_shifted.ResponseData(2,1,f)*exp(-1i*phase_shift);
end
Z_measured_gfl=Z_shifted;
Z_analytical_gfl=Z_wf;


% GFM
run('initialize_WT_default_gfm.m');
run('initialize_params_agg.m');
% cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
load('vals_pn_adm_WF_gfm_19062025.mat');

V_pn=vals{1};
I1_pn=vals{3};
sig=vals{5};

clear Y_c Z_c
for f=1:length(sig(1).freqs)
    Z_c(:,:,f)=[V_pn{1,f}, V_pn{2,f}]*([I1_pn{1,f}, I1_pn{2,f}]^-1);
end
Z_conv_exp=frd(-Z_c,sig(1).freqs);
Z_conv_exp.FrequencyUnit='Hz';

Z_shifted=Z_conv_exp;
phase_shift=-1.2;

for f=1:length(Z_shifted.Frequency)
    Z_shifted.ResponseData(1,2,f)=Z_shifted.ResponseData(1,2,f)*exp(1i*phase_shift);
    Z_shifted.ResponseData(2,1,f)=Z_shifted.ResponseData(2,1,f)*exp(-1i*phase_shift);
end
Z_measured_gfm=Z_shifted;
Z_analytical_gfm=Z_wf;
% 
% figure;
%%
Fs1=20;
f=figure;
f.Position=f.Position.*[00.2 0.2 2 2];
t=tiledlayout(size(Z_analytical_gfl.ResponseData,1)*2,size(Z_analytical_gfl.ResponseData,2),'TileSpacing',"tight",'Padding','compact');
mag=abs(Z_analytical_gfl.ResponseData);
mag_meas=abs(Z_measured_gfl.ResponseData);
P = atan2(real(Z_analytical_gfl.ResponseData),imag(Z_analytical_gfl.ResponseData));
phase=unwrap(P,[],3)*180/pi;
P2 = atan2(real(Z_measured_gfl.ResponseData),imag(Z_measured_gfl.ResponseData));
phase_meas=unwrap(P2,[],3)*180/pi;

mag2=abs(Z_analytical_gfm.ResponseData);
mag_meas2=abs(Z_measured_gfm.ResponseData);
P = atan2(real(Z_analytical_gfm.ResponseData),imag(Z_analytical_gfm.ResponseData));
phase2=unwrap(P,[],3)*180/pi;
P2 = atan2(real(Z_measured_gfm.ResponseData),imag(Z_measured_gfm.ResponseData));
phase_meas2=unwrap(P2,[],3)*180/pi;

titles={'$Z_{pp}$','$Z_{np}$';'$Z_{pn}$','$Z_{nn}$'};

for n_i=1:size(Z_analytical_gfl.ResponseData,1)
    for n_o=1:size(Z_analytical_gfl.ResponseData,2)
        num_m=tilenum(t,n_o*2-1,n_i);
        nexttile(num_m);
        m=loglog(Z_analytical_gfl.Frequency,squeeze(mag(n_o,n_i,:)),'Color',[0.0000 0.4470 0.7410]);
        grid on;hold on
        m2=loglog(Z_analytical_gfm.Frequency,squeeze(mag2(n_o,n_i,:)),'Color',[0.8500 0.3250 0.0980]);
        m.LineWidth=1;
        m2.LineWidth=1;
        mm=loglog(Z_measured_gfl.Frequency,squeeze(mag_meas(n_o,n_i,:)),'*','Color',[0.0000 0.4470 0.7410]);
        mm2=loglog(Z_measured_gfm.Frequency,squeeze(mag_meas2(n_o,n_i,:)),'*','Color',[0.8500 0.3250 0.0980]);
        
        set(gca,'FontSize',Fs1*0.7);
        if n_i==1
          ylabel('Mag (abs)','interpreter','latex','FontSize',Fs1);
        end
        set(gca,'XTickLabel',[]);
        if n_i==2
            % set(gca,'YTickLabel',[]);
        end
        if n_i==1 && n_o==1
            legend('GFL','GFM','interpreter','latex','Location','southwest','FontSize',Fs1);
        end

        title(titles{n_o,n_i},'Interpreter','latex','FontSize',Fs1);
        xlim([Z_measured_gfl.Frequency(1) Z_measured_gfl.Frequency(end)]);
        ylim([10^-5 10^3]);
    
        num_ph=tilenum(t,n_o*2,n_i);
        nexttile(num_ph);
        p=semilogx(Z_analytical_gfl.Frequency,squeeze(phase(n_o,n_i,:)),'Color',[0.0000 0.4470 0.7410]);
        grid on;hold on;
        p2=semilogx(Z_analytical_gfm.Frequency,squeeze(phase2(n_o,n_i,:)),'Color',[0.8500 0.3250 0.0980]);
        p.LineWidth=1;
        p2.LineWidth=1;
        pm=semilogx(Z_measured_gfl.Frequency,squeeze(phase_meas(n_o,n_i,:)),'*','Color',[0.0000 0.4470 0.7410]);
        pm2=semilogx(Z_measured_gfm.Frequency,squeeze(phase_meas2(n_o,n_i,:)),'*','Color',[0.8500 0.3250 0.0980]);
        
        xlim([Z_measured_gfl.Frequency(1) Z_measured_gfl.Frequency(end)]);
        ylim([-180 270]);
        % yticks([-180 -90 0 90 180]);
        set(gca,'FontSize',Fs1*0.7);
        if n_o==2
            xlabel('Frequency (Hz)','interpreter','latex','FontSize',Fs1);
        end
        if n_o==1
            set(gca,'XTickLabel',[]);
        end
        if n_i==1
            ylabel('Phase (deg)','interpreter','latex','FontSize',Fs1);
        end
        if n_i==2
            % set(gca,'YTickLabel',[]);
        end


    end
end

% exportgraphics(f,[path_for_save,'Z_measured','.pdf']);



%% Shift phase of impedance
% Z_shifted=Z_wf;
Z_shifted=Z_wt/(V_trans/V_coll)^2;

phase_shift=50*pi/180; %rad
% phase_shift=0.74; %rad
phase_shift=0.86; %rad

for f=1:length(Z_shifted.Frequency)
    Z_shifted.ResponseData(1,2,f)=Z_shifted.ResponseData(1,2,f)*exp(1i*phase_shift);
    Z_shifted.ResponseData(2,1,f)=Z_shifted.ResponseData(2,1,f)*exp(-1i*phase_shift);
end

%% Plot dd impedance with negative zone

% [Z_wf,Z_wt,cvtr,Z_wf_p1]=IM_WF(freq,Sb_WF,n_s,V_wt,P_pu,Q_pu,phi_wt,D_WW,cable_params,tf_params,Z_cable_g);

% Z=Z_wf_cable_g2;
% Z=Z_wf;
% Z=Z_C_C;
% Z=Z_line_G;
% Z=Z_wf_cable;
% Z=Z_wt;

out=1;
in=1;

f=figure;
% f=gcf;
f.Position=f.Position.*[0,0,1.5,0.8];
clear mag;
% mag=abs(Z_wf_cable_g.ResponseData);
mag=abs(Z2.ResponseData);
% mag2=abs(Z_eq_res_nl.ResponseData);
% mag3=abs(Z_eq_res.ResponseData);
% pass=real(Z_wf.ResponseData);
% pass=real(Z.ResponseData);
im_part=imag(Z2.ResponseData);

    %MIMO passivity
Z_dd=Z.ResponseData(1,1,:);
Z_qd=Z.ResponseData(1,2,:);
Z_dq=Z.ResponseData(2,1,:);
Z_qq=Z.ResponseData(2,2,:);

a=Z_dd+conj(Z_dd);
b=Z_qq+conj(Z_qq);
c=Z_dq+conj(Z_qd);

PD=Z;
PD.ResponseData=[a conj(c); c b];

for fr=1:length(PD.Frequency)
    PD_f=PD.ResponseData(:,:,fr);
    sigma=eig(PD_f);
    eig_f(:,fr)=sigma;
end

pass=min(eig_f,[],1);
% pass=real(Z(out,in).ResponseData);


fout=Z.Frequency;
p=loglog(fout,squeeze(mag(out,in,:)));
% p=semilogx(fout,squeeze(im_part(out,in,:)));
grid on;hold on;
% p2=loglog(fout,squeeze(mag2(out,in,:)));
% p3=loglog(fout,squeeze(mag3(out,in,:)));
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$\Re[Z_{pp}(s)] (\Omega)$','Interpreter','latex');
p.LineWidth=1.5;

[~,neg_ind]=find(pass<0.00);
ax=gca;
% ylim([1 ax.YLim(2)])
y_off=-0.02;
f_size=16;
ax.FontSize=f_size;
try
    for i=1:length(neg_ind)
        if i==1
            ind1=i;
            plot(fout(neg_ind(ind1)), 0,'k.');
            % text(fout(neg_ind(ind1)), ax.YLim(1)*1.2,[num2str(round(fout(neg_ind(ind1)),0)), ' Hz'],'HorizontalAlignment','right','VerticalAlignment','baseline','Interpreter','latex');
            text(fout(neg_ind(ind1)), y_off,[num2str(round(fout(neg_ind(ind1)),0)), ' Hz'],'HorizontalAlignment','left','VerticalAlignment','baseline','Interpreter','latex','FontSize',f_size);
            continue;
        else
            if neg_ind(i)-neg_ind(i-1)~=1
                ind1=i;
                plot(fout(neg_ind(ind1)), 0,'k.');
                % text(fout(neg_ind(ind1)), ax.YLim(1)*1.2,[num2str(round(fout(neg_ind(ind1)),0)), ' Hz'],'HorizontalAlignment','right','VerticalAlignment','baseline','Interpreter','latex');
                text(fout(neg_ind(ind1)),y_off ,[num2str(round(fout(neg_ind(ind1)),0)), ' Hz'],'HorizontalAlignment','right','VerticalAlignment','baseline','Interpreter','latex','FontSize',f_size);
%                 ind2=i;
            end
        end

        if i==length(neg_ind)
            ind2=i;
            text(fout(neg_ind(ind2)), y_off,[num2str(round(fout(neg_ind(ind2)),0)), ' Hz'],'HorizontalAlignment','right','VerticalAlignment','baseline','Interpreter','latex','FontSize',f_size);
            p=patch([fout(neg_ind(ind1)),fout(neg_ind(ind1)),fout(neg_ind(ind2)),fout(neg_ind(ind2)),fout(neg_ind(ind1))],[ax.YLim(1), ax.YLim(2), ax.YLim(2), ax.YLim(1), ax.YLim(1)],'r','FaceAlpha',0.3);
            plot(fout(neg_ind(ind2)), 0,'k.');
            % text(fout(neg_ind(ind2)), ax.YLim(1)*1.2,[num2str(round(fout(neg_ind(ind2)),0)), ' Hz'],'HorizontalAlignment','left','VerticalAlignment','baseline','Interpreter','latex');
           
        else
            if neg_ind(i+1)-neg_ind(i)==1
                continue;
            else
                ind2=i;
                text(fout(neg_ind(ind2)), y_off,[num2str(round(fout(neg_ind(ind2)),0)), ' Hz'],'HorizontalAlignment','right','VerticalAlignment','baseline','Interpreter','latex','FontSize',f_size);
                patch([fout(neg_ind(ind1)),fout(neg_ind(ind1)),fout(neg_ind(ind2)),fout(neg_ind(ind2)),fout(neg_ind(ind1))],[ax.YLim(1), ax.YLim(2), ax.YLim(2), ax.YLim(1), ax.YLim(1)],'r','FaceAlpha',0.3);
                plot(fout(neg_ind(ind2)), 0,'k.');
                % text(fout(neg_ind(ind2)), ax.YLim(1)*1.1,[num2str(round(fout(neg_ind(ind2)),0)), ' Hz'],'HorizontalAlignment','left','VerticalAlignment','baseline','Interpreter','latex');

            end
        end
    end  
end
[min_mag,min_ind]=min(mag(out,in,:));
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
p2=patch([fout(neg_ind(ind1)),fout(neg_ind(ind1)),fout(neg_ind(ind2)),fout(neg_ind(ind2)),fout(neg_ind(ind1))],[ax.YLim(1), ax.YLim(2), ax.YLim(2), ax.YLim(1), ax.YLim(1)],'r','FaceAlpha',0.3);
delete(p);
ylim([ax.YLim(1), ax.YLim(2)]);
% text(fout(min_ind)*1.06,min_mag*1.1,[num2str(round(fout(min_ind),0)), ' Hz'],'HorizontalAlignment','left','VerticalAlignment','baseline','Interpreter','latex');
% plot(fout(min_ind),min_mag,'k.');
% exportgraphics(f,[path_for_save,'Z_wf_passivity','.pdf']);



%% Plot WF and cable impedance with negative zone, vary cable distance

D_C=10:4:150;

f=figure;
t=tiledlayout(2,1);
t.TileSpacing='compact';
f.Position=f.Position.*[0,0,0.75*2,1.1];
f_size=16;

out=1;
in=1;
fout=Z_wf.Frequency;
p_vary=D_C;
color_plot = colormap(jet(length(p_vary)));
for i=1:length(p_vary)
    Z_line_C=IM_line(freq,r_c_eff*D_C(i),l_c_eff*D_C(i));    %Cable series imp
    Z_C_C=IM_C(freq,c_c_eff/2*D_C(i));                  %Cable shunt imp
    l_comp=(1/(c_c_eff/2*D_C(i)*2*pi*f1))/Z_b_G-(l_tf_wf+(l_c_eff*D_C(i)*2*pi*f1)/Z_b_G);
    L_comp=sqrt(2)*l_comp*(Z_b_G/(2*pi*f1));
    R_comp=0.1*l_comp*Z_b_G;
    Z_comp = IM_line(freq,R_comp,L_comp); 
    Z_line_OH=IM_line(freq,r_OH*Z_b_G*D_C(i)/50,l_OH*((Z_b_G/(2*pi*f1)))*D_C(i)/50);

     %Onshore case
    Z_eq=Z_line_G+Z_line_OH;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    Z_wf_OH_ser=Z_wf+Z_line_OH+Z_line_G;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;

    %Offshore case
    Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    % Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    % Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    % Z=Z_wf_cable_g;
    Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_wf_cable_g4=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;

    Z=Z_wf;
    Z2=Z_wf_OH;
    Z3=Z_wf_cable_g4;
    clear mag;
    clear pass;
        %MIMO passivity
  

    pass=mimo_passivity(Z);
    mag=abs(Z2.ResponseData);
    mag2=abs(Z3.ResponseData);
    % pass=real(Z(out,in).ResponseData);
    im_part=imag(Z2(out,in).ResponseData);
    im_part2=imag(Z3(out,in).ResponseData);
    figure(f);

    % nexttile(t,1);
    % p1=semilogx(fout,pass,'Color',color_plot(i,:));
    % grid on;hold on;
    % p1.LineWidth=1.5;

    nexttile(t,1);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;

    nexttile(t,2);
    p3=loglog(fout,squeeze(mag2(out,in,:)),'Color',color_plot(i,:));
    % p3=semilogx(fout,squeeze(im_part2(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p3.LineWidth=1.5;

end


% nexttile(t,1);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
% ax=gca;
% f_size=16;
% ax.FontSize=f_size;
% set(gca,'XTickLabel',[]);
% xlim([fout(1), fout(end)]);
% % ylim([-1200 200]);
% ylabel('$eig(Z_{wf}+Z_{wf}^H) (\Omega)$','Interpreter','latex');

nexttile(t,1);
ax=gca;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
line([neg_pass_freq neg_pass_freq],[ax.YLim(1) ax.YLim(end)],'Color','k','LineStyle','--');
xlim([fout(1), fout(end)]);
ylim([-50 320]);
% ylabel('$\Im\{Z_{on,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{on,pp}(j\omega)|$','Interpreter','latex');

nexttile(t,2);
ax=gca;
ax.FontSize=f_size;
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
line([neg_pass_freq neg_pass_freq],[ax.YLim(1) ax.YLim(end)],'Color','k','LineStyle','--');
xlim([fout(1), fout(end)]);
ylim([-2000 2000]);
% ylabel('$\Im\{Z_{off,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{off,pp}(j\omega)|$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');%,'Ticks',p_vary);
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'Distance (km)','Interpreter','latex');
% exportgraphics(f,[path_for_save,'vary_cable_gfl','.pdf']);
% exportgraphics(f,[path_for_save,'vary_cable_gfm','.pdf']);

%% Plot WF and cable impedance with negative zone, vary cable thickness

cable_data=xlsread('CableData.xlsx','HVAC');
cable_220=cable_data(5:8,:);

% [Z_wf,Z_wt,cvtr,Z_wf_p1]=IM_WF(freq,Sb_WF,n_s,V_wt,P_pu,Q_pu,phi_wt,D_WW,cable_params,tf_params,Z_cable_g);

%Default wind farm parameters
% Sb_WF=(12*1)*Sb_WT;
[Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;

D_C=50;

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,1];

out=1;
in=1;
fout=Z_wf.Frequency;
p_vary=cable_220;
color_plot = colormap(jet(size(p_vary,1)));
for i=1:size(p_vary,1)
    r_c=cable_220(i,4)*1e-3;      %Ohm/km
    l_c=cable_220(i,5)*1e-3;   %H/km
    c_c=cable_220(i,6)*1e-9;   %F/km
    
    n_cables=ceil(Sb_WF/(cable_220(i,8)*1e6));
    % n_cables=1;
    r_c_eff=r_c/n_cables;
    l_c_eff=l_c/n_cables;
    c_c_eff=c_c*n_cables;

    Z_line_C=IM_line(freq,r_c_eff*D_C,l_c_eff*D_C);    %Cable series imp
    Z_C_C=IM_C(freq,c_c_eff/2*D_C);                  %Cable shunt imp
    L_comp=sqrt(2)*((1/(c_c_eff/2*D_C*2*pi*f1))-(l_tf_wf*(Z_b_G)+(l_c_eff*D_C*2*pi*f1)))/(2*pi*f1);
    Z_comp = IM_line(freq,R_comp,L_comp); 

    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    Z_wf_line_g=Z_line_C+Z_wf+Z_line_G;
    Z_wf_cable_g2=Z_wf+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    Z_wf_cable_g3_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g4=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1)+Z_comp^-1)^-1;
    Z=Z_wf;
    Z2=Z_wf_cable_g3_ser;
    % Z=Z_wf_line_g;
    % Z=Z_wf;
    clear mag;
    clear pass;
        %MIMO passivity
    Z_dd=Z.ResponseData(1,1,:);
    Z_qd=Z.ResponseData(1,2,:);
    Z_dq=Z.ResponseData(2,1,:);
    Z_qq=Z.ResponseData(2,2,:);

    a=Z_dd+conj(Z_dd);
    b=Z_qq+conj(Z_qq);
    c=Z_dq+conj(Z_qd);

    PD=Z;
    PD.ResponseData=[a conj(c); c b];

    for fr=1:length(PD.Frequency)
        PD_f=PD.ResponseData(:,:,fr);
        sigma=eig(PD_f);
        eig_f(:,fr)=sigma;
    end
    % eig_A{i} = eig(cvtr.A);

    mag=abs(Z2.ResponseData);
    % pass=real(Z(out,in).ResponseData);
    pass=min(eig_f,[],1);
    im_part=imag(Z2(out,in).ResponseData);
    figure(f);

    % p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;
    if n_cables==4
        p2.LineStyle='--';
    end
    % [~,neg_ind]=find(pass(1,1,:)<0.00);
    [~,neg_ind]=find(pass<0.00);
    % mag_np=NaN(size(squeeze(mag(1,1,:))));
    % mag_np(neg_ind)=squeeze(mag(1,1,neg_ind));
    mag_np=NaN(size(squeeze(im_part(1,1,:))));
    mag_np(neg_ind)=squeeze(im_part(1,1,neg_ind));
    p_np=loglog(fout,mag_np,'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end
end

f_size=16;
ax=gca;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$\Im{Z_{pp}(s)}(\Omega)$','Interpreter','latex');



c=colorbar('FontSize',12,'TickLabelInterpreter','latex','Ticks',cable_220(:,2));
% c.Layout.Tile = 'east';
clim([cable_220(1,2) cable_220(end,2)]);
title(c,'Cable Size (mm2)','Interpreter','latex');
% exportgraphics(f,[path_for_save,'vary_cable','.pdf']);


%% Plot WF and cable impedance with negative zone, vary # turbines

Z=Z_wf;
% Sb_WF=[3:4:42]*Sb_WT;
% n_t=(3:4:42);

n_all=(Sb_WF/Sb_WT);
n_t=(n_all/2:4:n_all);
% n_t=200;

% Z_b_WF=(V_coll/sqrt(3))^2/Sb_WF;
% R_tf_wf=r_tf_wf*Z_b_WF;
% L_tf_wf=l_tf_wf*(Z_b_WF/(2*pi*f1));

f=figure;
t=tiledlayout(3,1);
t.TileSpacing='compact';
f.Position=f.Position.*[0,0,1.5,1.3];


out=1;
in=1;
fout=Z.Frequency;
p_vary=n_t;
p_ones=ones(1,length(p_vary));
params=[n_t.*p_ones;Sb_WT.*p_ones;D_WW.*p_ones];
% params=[Sb_WF.*p_ones;Sb_WT.*p_ones;n_s.*p_ones;D_WW.*p_ones];
color_plot = colormap(jet(length(p_vary)));
for i=1:length(p_vary)

    % tf_params=[R_tf_wt,L_tf_wt;R_tf_wf,L_tf_wf];
    [Z_wf,Z_conv,cvtr]=IM_WF_nt(freq,params(1,i),params(2,i),cvtr,params(3,i),cable_params,tf_params,V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    % Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    % Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    Z_wf_cable_g4=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1)+Z_comp^-1)^-1;
    % Z_wf_OH_ser=Z_wf+Z_line_OH+Z_line_G;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
    % Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
   

    Z=Z_wf;
    Z2=Z_wf_OH;
    Z3=Z_wf_cable_g4;
    % Z2=Z_wf;
    clear mag;
    clear pass;

    pass=mimo_passivity(Z);
    mag=abs(Z2.ResponseData);
    mag2=abs(Z3.ResponseData);
    im_part=imag(Z2(out,in).ResponseData);
    im_part2=imag(Z3(out,in).ResponseData);
    figure(f);

    nexttile(t,1);
    p1=semilogx(fout,pass,'Color',color_plot(i,:));
    grid on;hold on;
    p1.LineWidth=1.5;

    nexttile(t,2);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;

    nexttile(t,3);
    p3=loglog(fout,squeeze(mag2(out,in,:)),'Color',color_plot(i,:));
    % p3=semilogx(fout,squeeze(im_part2(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p3.LineWidth=1.5;
end


nexttile(t,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-1200 200]);
ylabel('Passivity','Interpreter','latex');

nexttile(t,2);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-1000 2000]);
% ylabel('$\Im\{Z_{on,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{on,pp}(j\omega)|$','Interpreter','latex');

nexttile(t,3);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
ax.FontSize=f_size;
% ylim([-1000 2000]);
xlim([fout(1), fout(end)]);
% ylabel('$\Im\{Z_{off,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{off,pp}(j\omega)|$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1)/n_all*100 p_vary(end)/n_all*100]);
title(c,'\% Turbines','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);
% exportgraphics(f,[path_for_save,'vary_turbines_gfl','.pdf']);
%exportgraphics(f,[path_for_save,'vary_turbines_gfm','.pdf']);

%% WF impedance with VSM control parameter variation, Active P ref
% Sb_WT=5e6;
% Sb_WF=(3*1)*5e6;
% n_s=1;
% V_wt=900;
P_pu=0.1:0.1:1;
Q_pu=0;
D_WW=1;


f=figure;
t=tiledlayout(3,1);
t.TileSpacing='compact';
f.Position=f.Position.*[0,0,1.5,1.3];
p_vary=P_pu;

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
for i=1:length(p_vary)
    cvtr.ss.pref0=p_vary(i);

    %Onshore case
    Z_eq=Z_line_G+Z_line_OH;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    % Z_wf_OH_ser=Z_wf+Z_line_OH+Z_line_G;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;

    %Offshore case
    Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    % Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    % Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    % Z=Z_wf_cable_g;
    Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_wf_cable_g4=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
   

    Z=Z_wf;
    Z2=Z_wf_OH;
    Z3=Z_wf_cable_g4;
    % Z2=Z_wf;
    clear mag;
    clear pass;
        %MIMO passivity
  

    pass=mimo_passivity(Z);
    mag=abs(Z2.ResponseData);
    mag2=abs(Z3.ResponseData);
    im_part=imag(Z2(out,in).ResponseData);
    im_part2=imag(Z3(out,in).ResponseData);
    figure(f);

    nexttile(t,1);
    p1=semilogx(fout,pass,'Color',color_plot(i,:));
    grid on;hold on;
    p1.LineWidth=1.5;

    nexttile(t,2);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;

    nexttile(t,3);
    p3=loglog(fout,squeeze(mag2(out,in,:)),'Color',color_plot(i,:));
    % p3=semilogx(fout,squeeze(im_part2(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p3.LineWidth=1.5;
    % [~,neg_ind]=find(pass(1,1,:)<0.00);
    % [~,neg_ind]=find(pass<0.00);
    % % mag_np=NaN(size(squeeze(mag(1,1,:))));
    % % mag_np(neg_ind)=squeeze(mag(1,1,neg_ind));
    % mag_np=NaN(size(squeeze(im_part(1,1,:))));
    % mag_np(neg_ind)=squeeze(im_part(1,1,neg_ind));
    % p_np=loglog(fout,mag_np,'Color','#515151');
    % try
    %     p_np.LineWidth=1.5;
    % end

    % nexttile(t2,3);
    % p3=semilogx(fout,squeeze(phase(1,1,:)),'Color',color_plot(i,:));
    % p3.LineWidth=1;
    % grid on;hold on;
    % yticks([-180 -90 0 90 180]);
end


nexttile(t,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-1200 200]);
ylabel('Passivity','Interpreter','latex');

nexttile(t,2);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-1000 2000]);
% ylabel('$\Im\{Z_{on,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{on,pp}(j\omega)|$','Interpreter','latex');

nexttile(t,3);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
% ylabel('$\Im\{Z_{off,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{off,pp}(j\omega)|$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');

% nexttile(t2,2);
% ax=gca;
% ax.FontSize=f_size;
% xlim([fout(1), fout(end)]);
% xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');
% 

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$P_{ref}$','interpreter','latex');
% exportgraphics(f,[path_for_save,'p_vary_gfl_zoom','.pdf']);
% exportgraphics(f,[path_for_save,'p_vary_gfm','.pdf']);

%% WF impedance with VSM control parameter variation, Reactive Q ref
% Sb_WT=5e6;
% Sb_WF=(3*1)*5e6;
% n_s=1;
% V_wt=900;
P_pu=0.5;
Q_pu=0:0.1:1;
% phi_wt=0;
% D_WW=1;
cvtr.ss.pref0=P_pu;

f=figure;
t=tiledlayout(3,1);
t.TileSpacing='compact';
% f=gcf;
f.Position=f.Position.*[0,0,1.5,1.3];

p_vary=Q_pu;

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
for i=1:length(p_vary)
    cvtr.ss.qref0=p_vary(i);
  
       %Onshore case
    Z_eq=Z_line_G+Z_line_OH;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    % Z_wf_OH_ser=Z_wf+Z_line_OH+Z_line_G;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;

    %Offshore case
    Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    % Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    % Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    % Z=Z_wf_cable_g;
    Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_wf_cable_g4=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
   

    Z=Z_wf;
    Z2=Z_wf_OH;
    Z3=Z_wf_cable_g4;
    % Z2=Z_wf;
    clear mag;
    clear pass;
        %MIMO passivity
  

    pass=mimo_passivity(Z);
    mag=abs(Z2.ResponseData);
    mag2=abs(Z3.ResponseData);
    im_part=imag(Z2(out,in).ResponseData);
    im_part2=imag(Z3(out,in).ResponseData);
    figure(f);

    nexttile(t,1);
    p1=semilogx(fout,pass,'Color',color_plot(i,:));
    grid on;hold on;
    p1.LineWidth=1.5;

    nexttile(t,2);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;

    nexttile(t,3);
    p3=loglog(fout,squeeze(mag2(out,in,:)),'Color',color_plot(i,:));
    % p3=semilogx(fout,squeeze(im_part2(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p3.LineWidth=1.5;
    % [~,neg_ind]=find(pass(1,1,:)<0.00);
    % [~,neg_ind]=find(pass<0.00);
    % % mag_np=NaN(size(squeeze(mag(1,1,:))));
    % % mag_np(neg_ind)=squeeze(mag(1,1,neg_ind));
    % mag_np=NaN(size(squeeze(im_part(1,1,:))));
    % mag_np(neg_ind)=squeeze(im_part(1,1,neg_ind));
    % p_np=loglog(fout,mag_np,'Color','#515151');
    % try
    %     p_np.LineWidth=1.5;
    % end

    % nexttile(t2,3);
    % p3=semilogx(fout,squeeze(phase(1,1,:)),'Color',color_plot(i,:));
    % p3.LineWidth=1;
    % grid on;hold on;
    % yticks([-180 -90 0 90 180]);
end


nexttile(t,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-1200 200]);
ylabel('Passivity','Interpreter','latex');

nexttile(t,2);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-1000 2000]);
% ylabel('$\Im\{Z_{on,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{on,pp}(j\omega)|$','Interpreter','latex');

nexttile(t,3);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
% ylabel('$\Im\{Z_{off,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{off,pp}(j\omega)|$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');

% nexttile(t2,2);
% ax=gca;
% ax.FontSize=f_size;
% xlim([fout(1), fout(end)]);
% xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');
% 

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$Q_{ref}$','interpreter','latex');
% exportgraphics(f,[path_for_save,'q_vary_gfl_zoom','.pdf']);
% exportgraphics(f,[path_for_save,'q_vary_gfm','.pdf']);

%% WF impedance with VSM control parameter variation, grid impedance
% Sb_WT=5e6;
% Sb_WF=(12*1)*5e6;
% n_s=1;
% V_wt=900;
% P_pu=1;
% Q_pu=0;
% phi_wt=0;
% D_WW=1;
% %Rg
% rg=[0.001, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03,0.04,0.07, 0.1];
% rg=linspace(0.001,0.1,9);
% lg=0.1;
% cvtr.parameters.lg=lg;
% Lg
rg=0.01;
cvtr.parameters.rg=rg;
% lg=[0.001, 0.005, 0.01, 0.015, 0.02,0.05, 0.1,0.2,0.3,0.4];
lg=linspace(0.01, 1,12);

p_vary=lg;

f=figure;
t=tiledlayout(2,1);
t.TileSpacing='compact';
f.Position=f.Position.*[0,0,0.75*2,1.1];

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));

for i=1:length(p_vary)
    if p_vary==rg
        cvtr.parameters.rg=p_vary(i);    
        Rg=p_vary(i)*Z_b_G;
        Lg=lg*((Z_b_G/(2*pi*f1)));
    elseif p_vary==lg
        cvtr.parameters.lg=p_vary(i);    
        Rg=rg*Z_b_G;
        Lg=p_vary(i)*((Z_b_G/(2*pi*f1)));
    end

    Z_line_G=IM_line(freq,Rg,Lg);    %Grid equiv.

       %Onshore case
    Z_eq=Z_line_G+Z_line_OH;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;
    % Z_wf_OH_ser=Z_wf+Z_line_OH+Z_line_G;
    Z_wf_OH=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;

    %Offshore case
    Z_eq=(((Z_comp^-1+Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_f1=Z_eq.ResponseData(:,:,1);
    r_g_agg=(r_tf_wt+r_tf_wf+real(Z_f1(1,1)/Z_b_G)/3);
    l_g_agg=(l_tf_wt+l_tf_wf+imag(Z_f1(1,1)/Z_b_G)/3);
    cvtr_agg=initialize_WT_default_agg(cvtr,Sb_WF,r_g_agg,l_g_agg);
    [Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    % Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    % Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    % Z=Z_wf_cable_g;
    Z_wf_cable_g4_ser=Z_wf+(((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;
    Z_wf_cable_g4=(Z_wf^-1+((Z_C_C^-1+Z_line_G^-1+Z_comp^-1)^-1+Z_line_C)^-1+Z_C_C^-1+Z_comp^-1)^-1;

    
    Z=Z_wf;
    Z2=Z_wf_OH;
    Z3=Z_wf_cable_g4;
    clear mag;
    clear pass;
        %MIMO passivity
  

    pass=mimo_passivity(Z);
    mag=abs(Z2.ResponseData);
    mag2=abs(Z3.ResponseData);
    im_part=imag(Z2(out,in).ResponseData);
    im_part2=imag(Z3(out,in).ResponseData);
    figure(f);

    % nexttile(t,1);
    % p1=semilogx(fout,pass,'Color',color_plot(i,:));
    % grid on;hold on;
    % p1.LineWidth=1.5;

    nexttile(t,1);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(pass),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;

    nexttile(t,2);
    p3=loglog(fout,squeeze(mag2(out,in,:)),'Color',color_plot(i,:));
    % p3=semilogx(fout,squeeze(im_part2(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p3.LineWidth=1.5;

end


% nexttile(t,1);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
% ax=gca;
% f_size=16;
% ax.FontSize=f_size;
% set(gca,'XTickLabel',[]);
% xlim([fout(1), fout(end)]);
% % ylim([-1200 200]);
% ylabel('$eig(Z_{wf}+Z_{wf}^H) (\Omega)$','Interpreter','latex');

nexttile(t,1);
ax=gca;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
line([neg_pass_freq neg_pass_freq],[ax.YLim(1) ax.YLim(end)],'Color','k','LineStyle','--');
xlim([fout(1), fout(end)]);
% ylim([-50 320]);
% ylabel('$\Im\{Z_{on,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{on,pp}(j\omega)|$','Interpreter','latex');

nexttile(t,2);
ax=gca;
ax.FontSize=f_size;
% line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
line([neg_pass_freq neg_pass_freq],[ax.YLim(1) ax.YLim(end)],'Color','k','LineStyle','--');
xlim([fout(1), fout(end)]);
% ylim([-2000 2000]);
% ylabel('$\Im\{Z_{off,pp}(j\omega)\}$','Interpreter','latex');
ylabel('$|Z_{off,pp}(j\omega)|$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$L_g (pu)$','interpreter','latex');
% exportgraphics(f,[path_for_save,'vary_lg_gfl','.pdf']);
% exportgraphics(f,[path_for_save,'vary_lg_gfm','.pdf']);

%% Any control parameter

% run('initialize_WT_default.m');
clear p_vary
% cvtr.parameters.model_selector=3;
cvtr_nom=cvtr;
% kad=[0, 0.1, 0.5, 0.8, 1, 1.1, 2, 4];
%Active damping
kad=linspace(0, 4, 9);
Omegaad=linspace(1e-4,10,15)*2*pi;
% Tfad=linspace(1e-1,1e1,16);

%Feedforward
kffv=[0 1];

% Tfpll=[0.0002, .002, .02, .2];

%Virtual impedance
rs=linspace(0,.5,6);
ls=linspace(0.1,1,10);

rg=linspace(0,0.4,10)/0.01;
lg=linspace(0.01,0.8,10)/0.1;

%Droop and VSM parameters
ta = linspace(1, 100,10)/50;
kd = linspace(1,100,10)/40;
kdrpOmega = linspace(1,100,10)/20;
kdrpq = linspace(0.005,0.15,10)/0.05;
Omegadf = linspace(10,100,9)/50;
Omegaqf = linspace(10,100,9)/100;

%Filter parameters
cf = linspace(1e-5,1e-2,20);
lf = linspace(0.001,0.4,20);

p_vary.name='ls';
p_vary.val=eval(p_vary.name);

f=figure;

f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,1];


out=1;
in=1;
fout=Z_wt.Frequency;

color_plot = colormap(jet(length(p_vary.val)));
% color_plot = flip(colormap(jet(length(p_vary))));
nom=getfield(cvtr_nom,'parameters', p_vary.name);
for i=1:length(p_vary.val) 
    cvtr=setfield(cvtr,'parameters', p_vary.name,p_vary.val(i)*nom); 
    [Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
       Z_wf=Z_wf*(V_trans/V_coll)^2;

    Z_line_G=IM_line(freq,Rg,Lg);    %Grid equiv.
    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    Z=Z_wf;
    Z2=Z_wf_cable_g3;
    % Z=Z_wt;
    clear mag pass phase
    mag=abs(Z2(out,in).ResponseData);
    phase=angle(Z2(out,in).ResponseData)*180/pi;
    % pass=real(Z(out,in).ResponseData);
    pass=mimo_passivity(Z);
    im_part=imag(Z2(out,in).ResponseData);
     
    figure(f);
    % p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(pass(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;

    eig_A=eig(cvtr.A);
    if any(real(eig_A)>0)
        p1.LineStyle='--';
    end
    [~,neg_ind]=find(pass<0.00);
    % mag_np=NaN(size(squeeze(mag(1,1,:))));
    % mag_np(neg_ind)=squeeze(mag(1,1,neg_ind));
    mag_np=NaN(size(squeeze(im_part(1,1,:))));
    mag_np(neg_ind)=squeeze(im_part(1,1,neg_ind));
    % p_np=loglog(fout(neg_ind),squeeze(mag(out,in,neg_ind)),'Color','#515151');
    p_np=loglog(fout,mag_np,'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end
end




figure(f);
% nexttile(t2,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
% ylim([-200 300]);
ylabel('$\min\{\lambda_{1,2}[Z(s)+Z^H(s)]\} (pu)$','Interpreter','latex');
xlabel('Frequency (Hz)','interpreter','latex');

c=colorbar('FontSize',12,'TickLabelInterpreter','latex','Ticks',p_vary.val);
% c.Layout.Tile = 'east';
clim([p_vary.val(1) p_vary.val(end)]);
title(c,p_vary.name,'interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);


%% WF impedance with VSM control parameter variation, other control variables

% %Rg
% rg=[0.001, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03,0.04,0.07, 0.1];
% lg=0.01;
% cvtr.parameters.lg=lg;
%Lg

rs=[.015,.02,.03,0.04,0.05,0.06,0.1,0.2];
% ls=[0.01,0.05,0.1,0.2,0.3,0.4];

p_vary=ls;

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,0.8*2];
t2=tiledlayout(2,1);
t2.TileSpacing='tight';
nexttile(t2);

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
for i=1:length(p_vary)
    cvtr.parameters.rs=p_vary(i);    

    [Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    Z_line_G=IM_line(freq,Rg,Lg);    %Grid equiv.
    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z=Z_wf_cable_g;
    Z=Z_conv;
    clear mag;
    clear pass;
    mag=abs(Z.ResponseData);
    pass=real(Z(out,in).ResponseData);
    figure(f);
    nexttile(t2,1);
    p=semilogx(fout,squeeze(pass(1,1,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p.LineWidth=1.5;
    

    nexttile(t2,2);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;
    [~,neg_ind]=find(pass(out,in,:)<0.00);
    mag_np=NaN(size(squeeze(mag(out,in,:))));
    mag_np(neg_ind)=squeeze(mag(out,in,neg_ind));
    % p_np=loglog(fout(neg_ind),squeeze(mag(out,in,neg_ind)),'Color','#515151');
    p_np=loglog(fout,mag_np,'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end
end

nexttile(t2,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
ylim([-200 300]);
ylabel('$\Re\{Z_{pp}(s)\} (\Omega)$','Interpreter','latex');


nexttile(t2,2);
ax=gca;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');


c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$Param$','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);


%% WF impedance with VSM control parameter variation, filter impedance
Sb_WT=5e6;
Sb_WF=(3*1)*5e6;
n_s=1;
V_wt=900;
P_pu=1;
Q_pu=0;
phi_wt=0;
D_WW=1;

% %Change resistance
% rf=[0.001, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03,0.04,0.07, 0.1];
% lf=0.1;
% cf=7.4e-4;
% cvtr.parameters.lf=lf;
% cvtr.parameters.cf=cf;

% %Change inductance
% rf=.1;
% lf=[0.001, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03,0.04,0.07, 0.1, 0.2];
% cf=7.4e-4;
% cvtr.parameters.rf=rf;
% cvtr.parameters.cf=cf;

%Change capacitance
rf=.01;
lf=.1;
cf=logspace(-5,-2,10);
cvtr.parameters.rf=rf;
cvtr.parameters.lf=lf;

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,0.8*2];
t2=tiledlayout(2,1);
t2.TileSpacing='tight';
nexttile(t2);

p_vary=cf;

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
for i=1:length(p_vary)
    if p_vary==rf
        cvtr.parameters.rf=p_vary(i); 
    elseif p_vary==lf
        cvtr.parameters.lf=p_vary(i); 
    elseif p_vary==cf
        cvtr.parameters.cf=p_vary(i); 
    end


    [Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    Z_line_G=IM_line(freq,Rg,Lg);    %Grid equiv.
    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z=Z_wf_cable_g;
    Z=Z_conv;
    clear mag;
    clear pass;
    mag=abs(Z.ResponseData);
    pass=real(Z(out,in).ResponseData);
    figure(f);
    nexttile(t2,1);
    p=semilogx(fout,squeeze(pass(1,1,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p.LineWidth=1.5;
    

    nexttile(t2,2);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;
    [~,neg_ind]=find(pass(out,in,:)<0.00);
    p_np=loglog(fout(neg_ind),squeeze(mag(out,in,neg_ind)),'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end
end

nexttile(t2,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
ylim([-200 300]);
ylabel('$\Re\{Z_{pp}(s)\} (\Omega)$','Interpreter','latex');


nexttile(t2,2);
ax=gca;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');


c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$Param$','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);

%% WF impedance with VSM control parameter variation, any control parameter
run('initialize_WT_default.m');

% kad=[0, 0.1, 0.5, 0.8, 1, 1.1, 2, 4];
kad=linspace(0, 4, 9);
Omegaad=[0.001,0.01,0.1,10, 20, 50, 100, 1000]*2*pi;
kffv=[0, 0.01, 0.1, 0.2, 0.5, 1, 2, 4, 10];
Tfpll=[0.0002, .002, .02, .2];

p_vary.name='kad';
p_vary.val=eval(p_vary.name);

f=figure;

% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,0.8*2];
t2=tiledlayout(2,1);
t2.TileSpacing='tight';
nexttile(t2);

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary.val)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
for i=1:length(p_vary.val)
    
    cvtr=setfield(cvtr,'parameters', p_vary.name,p_vary.val(i)); 

    [Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    % Z=Z_wf_cable_g;
    Z=Z_conv;
    clear mag;
    clear pass;

    %MIMO passivity
    Z_dd=Z.ResponseData(1,1,:);
    Z_qd=Z.ResponseData(1,2,:);
    Z_dq=Z.ResponseData(2,1,:);
    Z_qq=Z.ResponseData(2,2,:);

    a=Z_dd+conj(Z_dd);
    b=Z_qq+conj(Z_qq);
    c=Z_dq+conj(Z_qd);

    PD=Z;
    PD.ResponseData=[a conj(c); c b];

    for fr=1:length(PD.Frequency)
        PD_f=PD.ResponseData(:,:,fr);
        % PD_f=Z.ResponseData(:,:,f)+Z.ResponseData(:,:,f)';
        % PD.ResponseData(:,:,f)=Z.ResponseData(:,:,f)+Z.ResponseData(:,:,f)';
        % if freq==675
        %     disp('');
        % end
        sigma=eig(PD_f);
        eig_f(:,fr)=sigma;
    end
    eig_A{i} = eig(cvtr.A);

    mag=abs(Z.ResponseData);
    % pass=real(Z(out,in).ResponseData);
    pass=min(eig_f,[],1);
    figure(f);
    nexttile(t2,1);
    p=semilogx(fout,squeeze(pass),'Color',color_plot(i,:));
    grid on;hold on;
    p.LineWidth=1.5;

    if any(real(eig_A{i})>0)
        p.LineStyle='--';
    end
    

    nexttile(t2,2);
    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;
    [~,neg_ind]=find(pass<0.00);
    mag_np=NaN(size(squeeze(mag(out,in,:))));
    mag_np(neg_ind)=squeeze(mag(out,in,neg_ind));
    % p_np=loglog(fout(neg_ind),squeeze(mag(out,in,neg_ind)),'Color','#515151');
    p_np=loglog(fout,mag_np,'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end
    
    if any(real(eig_A{i})>0)
        p2.LineStyle='--';
        p_np.LineStyle='--';
    end
end


figure(f);
nexttile(t2,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-200 300]);
ylabel('$\Re\{Z_{pp}(s)\} (\Omega)$','Interpreter','latex');


nexttile(t2,2);
ax=gca;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');


c=colorbar('FontSize',12,'TickLabelInterpreter','latex','Ticks',p_vary.val);
c.Layout.Tile = 'east';
clim([p_vary.val(1) p_vary.val(end)]);
title(c,'$Param$','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);


%Plot eigenvalues
f_eig=figure;
for i=1:length(p_vary.val)
    scatter(real(eig_A{i}),imag(eig_A{i}),'or','MarkerEdgeColor',color_plot(i,:));
    xlabel('Real','interpreter','latex');ylabel('Imag');
    % xlim([-500,0])
    grid on; hold on
end


%% Plot FFT of voltage/current
V=out.V_cable;
I=out.I_cable;

% V=out.V_g;
% I=out.I_g;

% % 
% V=out.V_conv;
% I=out.I_conv;
% % 
% V=out.V_dq
% I=out.I_dq;

Fs=1/(V.time(2)-V.time(1));
N=length(V.time);
freqs=round(Fs/N*(0:N-1),4);

V=fft(V.signals.values)/(N/2);
I=fft(I.signals.values)/(N/2);

figure;
% plot(freqs,abs(V_g),"LineWidth",1);
% hold on;
% plot(freqs,abs(V),"LineWidth",1);hold on
plot(freqs,abs(I),"LineWidth",1);hold on
xlim([0 2000]);