clear;
path_for_save='Figures/';
set(0,'defaulttextinterpreter','latex');
% run('initialize_WT_default_GFM.m');
run('initialize_WT_default_GFL_onshore_SSO.m');
% run('initialize_WT_default_SSO_GFL.m');
% run('initialize_WT_default_220125.m');

freq=logspace(-1,5,800);

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
% V_coll=33e3;
% V_trans=V_coll;

D_WW=1;

% Cable and grid impedance
Z_b_G=(V_trans/sqrt(3/2))^2/Sb_WF;
% Z_b_G=(V_trans)^2/(Sb_WT);
D_C=90;    %km
%Values for 220kV, 360MW cable
r_c=27e-3;      %Ohm/km
l_c=0.386e-3;   %H/km
% c_c=0.177e-6;
c_c=0.177e-6;   %F/km
% Rg=0.86;
% Lg=99.03e-3;
% Rg=0.001/Z_b_G;
% Lg=0.01/((Z_b_G/(2*pi*f1)));
% Rg=0.001*Z_b_G;
% Lg=0.01*((Z_b_G/(2*pi*f1)));
%Overhead line
r_OH=0.02;
l_OH=0.2;
R_OH=r_OH*Z_b_G;
L_OH=l_OH*((Z_b_G/(2*pi*f1)));


Rg=0.01*Z_b_G;
Lg=0.60*((Z_b_G/(2*pi*f1)));



% cvtr.parameters.rg=Rg/(V_wt^2/Sb_WT);
% cvtr.parameters.lg=Lg/(V_wt^2/Sb_WT/(2*pi*f1));

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
Z_line_OH=IM_line(freq,R_OH,L_OH);

% %Convert to 33kV
% Z_line_C=Z_line_C*(33e3/220e3)^2;
% Z_C_C=Z_C_C*(33e3/220e3)^2;
% Z_line_G=Z_line_G*(33e3/220e3)^2;

% Z_cable_g=(((Z_line_G^-1+Z_C_C^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;

% Z_test=(Z_line_C^-1+Z_C_C^-1)^-1;
% Z_test_ser=Z_C_C+Z_line_C;
% Z_test2=Z_test+Z_line_C;

% WF impedance

%Impedance params
r_ww=0.098;
l_ww=0.36e-3;
c_ww=0.23e-6;
r_sw=0.041;
l_sw=0.31e-3;
c_sw=0.34e-6;
cable_params=[r_ww,l_ww,c_ww;r_sw,l_sw,c_sw];

Z_b_WT=(V_coll/sqrt(3/2))^2/Sb_WT;
Z_b_WF=(V_coll/sqrt(3/2))^2/Sb_WF;

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


% Z_wf=Z_wf*(V_trans/V_coll)^2;
% run('initialize_adm_WT.m');
% run('initialize_WT.m');

r_g_agg=r_tf_wt+r_tf_wt;%+r_OH+Rg/Z_b_G;
l_g_agg=l_tf_wt+l_tf_wt;%+l_OH+Lg/(Z_b_G/(2*pi*f1));

run('initialize_WT_default_agg.m');

[Z_wf,Z_conv,cvtr_agg,Z_wt]=IM_WF_agg(freq,Sb_WF,cvtr_agg,[r_tf_wt,l_tf_wt;r_tf_wf,l_tf_wf],V_coll);
Z_trans_WF=IM_line(freq,r_tf_wf*Z_b_WF,l_tf_wf*(Z_b_WF/(2*pi*cvtr.pu.fb))); %At 33k

Z_wf=Z_wf*(V_trans/V_coll)^2;
Z_wt=Z_wt*(V_trans/V_coll)^2;
% WF and Cable

Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_cable_g=Z_wf_cable+Z_line_G;
Z_wf_cable_g2=Z_wf+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_cable2=Z_wf+((Z_C_C^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
% Z_wf_cable_g=(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1+Z_wf;

% Z_wf_cable_parall=(Z_wf^-1+Z_cable_g^-1)^-1;
Z_wf_line_g=(Z_wf^-1+(Z_line_OH+Z_line_G)^-1)^-1;
Z_wf_line_g_ser=Z_wf+Z_line_OH+Z_line_G;
% % % Impedance from a WT
% Z_wf_wt_p=(Z_wf_p1^-1+Z_wt^-1)^-1;
% Z_wf_wt_s=Z_wf_p1+Z_wt;


% Z_filt_L=IM_line(freq,cvtr.parameters.Rf,cvtr.parameters.Lf);    %Cable series imp
% Z_filt_C=IM_C(freq,cvtr.parameters.Cf);   
% Z_test=(Z_filt_L^-1+Z_filt_C^-1)^-1*cvtr.pu.Zb*(V_coll/cvtr.parameters.Vnr)^2;
% Z_trans_WT=IM_line(freq,R_tf_wt,L_tf_wt);
% Z_test2=Z_test+Z_trans_WT;

% %% Plot impedance
opts_bd=bodeoptions;
opts_bd.FreqUnits='Hz';
opts_bd.Grid='on';
opts_bd.MagUnits='abs';
opts_bd.MagScale='log';
opts_bd.IOGrouping='none';
% figure;
% Z_wt=Y_wt^-1;
% 
% figure;
% % bodeplot(Z_wf,opts_bd);
% bodeplot(Z_wf_cable_g2,opts_bd);
% % bodeplot(Z_wf_cable_g(1,1),opts_bd);
% % bodeplot(Z_wf_cable_parall(1,1),opts_bd);
% % bodeplot(Z1,opts_bd);hold on

% % bodeplot(Z2,opts_bd,'*');
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


%% Plot pass and imag part of Z
% Z=Z_wt;
% Z=Z_wf;
Z=Z_wf;
% Z=Z_wf_cable_g2;
Z2=Z_wf_line_g_ser;
% Z2=Z_wf;
% Z2=Z_wf_cable_g3;
% Z2=Z_wf;
% Z1=Z2;

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

out=1;
in=1;
fout=Z.Frequency;
nexttile(t);
p=semilogx(fout,pass);
grid on;
p.LineWidth=1.5;
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;
ylabel('Passivity');

nexttile(t);
p=semilogx(fout,squeeze(im_part(out,in,:)));
grid on
% p2=loglog(fout,squeeze(mag(out,in,:)));
% p3=loglog(fout,squeeze(mag3(out,in,:)));
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$\Re[Z_{pp}(s)] (\Omega)$','Interpreter','latex');
p.LineWidth=1.5;

ax=gca;
% ylim([1 ax.YLim(2)])
f_size=16;
ax.FontSize=f_size;


line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');



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

%% MIMO Passivity
Z=Z_conv;
% Z=Z_wf_cable_g;
% Z=Z_wf_cable;
% Z=Z2;

Z_dd=Z.ResponseData(1,1,:);
Z_qd=Z.ResponseData(1,2,:);
Z_dq=Z.ResponseData(2,1,:);
Z_qq=Z.ResponseData(2,2,:);

a=Z_dd+conj(Z_dd);
b=Z_qq+conj(Z_qq);
c=Z_dq+conj(Z_qd);

PD=Z;
PD.ResponseData=[a conj(c); c b];

for f=1:length(PD.Frequency)
    PD_f=PD.ResponseData(:,:,f);
    % PD_f=Z.ResponseData(:,:,f)+Z.ResponseData(:,:,f)';
    % PD.ResponseData(:,:,f)=Z.ResponseData(:,:,f)+Z.ResponseData(:,:,f)';
    if f==675
        disp('');
    end
    sigma=eig(PD_f);
    eig_f(:,f)=sigma;
end

fout=Z.Frequency;

figure;
t2=tiledlayout(2,1);
t2.TileSpacing='compact';
clear mag;
mag=eig_f;

% mag=1./real(Z.ResponseData);
for out=1:size(mag,1)
        nexttile(t2);
        p=semilogx(fout,squeeze(mag(out,:)));
    %     p=loglog(fout,squeeze(mag(out,1,:)));
        grid on;hold on;
        xlabel('Frequency (Hz)');ylabel([Z.OutputName{out}, ', Mag (abs)']);
    %     yticks(floor(min(phase(out,1,:))/90)*90:90:ceil(max(phase(out,1,:))/90)*90);
        
        if out==1
            title('Eig 1');%,'Interpreter','latex');
        else
            title('Eig 2');
        end
    
        [~,neg_ind]=find(mag(out,:)<0.00);
        ax=gca;
        for i=1:length(neg_ind)
            if i==1
                ind1=i;
                continue;
            else
                if neg_ind(i)-neg_ind(i-1)~=1
                    ind1=i;
    %                 ind2=i;
                end
            end
    
            if i==length(neg_ind)
                ind2=i;
                patch([fout(neg_ind(ind1)),fout(neg_ind(ind1)),fout(neg_ind(ind2)),fout(neg_ind(ind2)),fout(neg_ind(ind1))],[ax.YLim(1), ax.YLim(2), ax.YLim(2), ax.YLim(1), ax.YLim(1)],'r','FaceAlpha',0.3);
            else
                if neg_ind(i+1)-neg_ind(i)==1
                    continue;
                else
                    ind2=i;
                    patch([fout(neg_ind(ind1)),fout(neg_ind(ind1)),fout(neg_ind(ind2)),fout(neg_ind(ind2)),fout(neg_ind(ind1))],[ax.YLim(1), ax.YLim(2), ax.YLim(2), ax.YLim(1), ax.YLim(1)],'r','FaceAlpha',0.3);
                end
            end
        end
end


%% Plot WF and cable impedance with negative zone, vary cable distance

% [Z_wf,Z_wt,cvtr,Z_wf_p1]=IM_WF(freq,Sb_WF,n_s,V_wt,P_pu,Q_pu,phi_wt,D_WW,cable_params,tf_params,Z_cable_g);

%Default wind farm parameters
% Sb_WF=(12*1)*Sb_WT;
[Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
Z_wf=Z_wf*(V_trans/V_coll)^2;

D_C=10:4:100;

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,1];

out=1;
in=1;
fout=Z_wf.Frequency;
p_vary=D_C;
color_plot = colormap(jet(length(p_vary)));
for i=1:length(p_vary)
    Z_line_C=IM_line(freq,r_c_eff*D_C(i),l_c_eff*D_C(i));    %Cable series imp
    Z_C_C=IM_C(freq,c_c_eff/2*D_C(i));                  %Cable shunt imp
    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    Z_wf_line_g=Z_line_C+Z_wf+Z_line_G;
    Z_wf_cable_g2=Z_wf+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    Z=Z_wf;
    Z2=Z_wf_cable_g3;
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
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');


c=colorbar('FontSize',12,'TickLabelInterpreter','latex');%,'Ticks',p_vary);
% c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'Cable Distance (km)','Interpreter','latex');
% exportgraphics(f,[path_for_save,'vary_cable','.pdf']);

%% Plot WF and cable impedance with negative zone, vary # turbines

Z=Z_wf;
% Sb_WF=[3:4:42]*Sb_WT;
% n_t=(3:4:42);

n_t=(0:4:(Sb_WF/Sb_WT));

% Z_b_WF=(V_coll/sqrt(3))^2/Sb_WF;
% R_tf_wf=r_tf_wf*Z_b_WF;
% L_tf_wf=l_tf_wf*(Z_b_WF/(2*pi*f1));

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,1];

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

    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    Z=Z_wf;
    Z2=Z_wf_cable_g3;
       clear mag;
    clear pass;
        %MIMO passivity
    pass=mimo_passivity(Z);

    mag=abs(Z2.ResponseData);
    % pass=real(Z(out,in).ResponseData);
    im_part=imag(Z(out,in).ResponseData);
    figure(f);

    % p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;
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
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');



c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
% c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'\# Turbines','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);
exportgraphics(f,[path_for_save,'vary_turbines','.pdf']);

%% WF impedance with VSM control parameter variation, Active P ref
% Sb_WT=5e6;
% Sb_WF=(3*1)*5e6;
% n_s=1;
% V_wt=900;
P_pu=0.1:0.1:1;
Q_pu=0;
D_WW=1;

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,0.8*2];
% t2=tiledlayout(3,1);
% t2.TileSpacing='tight';
% nexttile(t2);

p_vary=P_pu;

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
for i=1:length(p_vary)
    cvtr.ss.pref0=p_vary(i);
    [Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    % Z=Z_wf_cable_g;
    Z=Z_wf;
    Z2=Z_wf_cable_g3;
    % Z2=Z_wf;
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
    im_part=imag(Z(out,in).ResponseData);
    figure(f);

    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;
    % [~,neg_ind]=find(pass(1,1,:)<0.00);
    [~,neg_ind]=find(pass<0.00);
    mag_np=NaN(size(squeeze(mag(1,1,:))));
    mag_np(neg_ind)=squeeze(mag(1,1,neg_ind));
    p_np=loglog(fout,mag_np,'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end

    % nexttile(t2,3);
    % p3=semilogx(fout,squeeze(phase(1,1,:)),'Color',color_plot(i,:));
    % p3.LineWidth=1;
    % grid on;hold on;
    % yticks([-180 -90 0 90 180]);
end


% nexttile(t2,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-200 300]);
ylabel('$\Re\{Z_{pp}(s)\} (\Omega)$','Interpreter','latex');


% nexttile(t2,2);
% ax=gca;
% ax.FontSize=f_size;
% xlim([fout(1), fout(end)]);
% xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');
% 

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
% c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$P_{ref}$','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);

%% WF impedance with VSM control parameter variation, Reactive Q ref
% Sb_WT=5e6;
% Sb_WF=(3*1)*5e6;
% n_s=1;
% V_wt=900;
P_pu=0.5;
Q_pu=0:0.1:1;
% phi_wt=0;
% D_WW=1;

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,0.8*2];
% t2=tiledlayout(3,1);
% t2.TileSpacing='tight';
% nexttile(t2);

p_vary=Q_pu;

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
for i=1:length(p_vary)
    cvtr.ss.qref0=p_vary(i);
    [Z_wf,Z_conv,cvtr]=IM_WF2(freq,Sb_WF,Sb_WT,n_s,cvtr,D_WW,cable_params,tf_params,V_coll);
    Z_wf=Z_wf*(V_trans/V_coll)^2;

    Z_wf_cable=(((Z_C_C^-1+Z_wf^-1)^-1+Z_line_C)^-1+Z_C_C^-1)^-1;
    Z_wf_cable_g=Z_wf_cable+Z_line_G;
    Z_wf_cable_g3=(Z_wf^-1+(((Z_C_C^-1+Z_line_G^-1)^-1+Z_line_C)^-1+Z_C_C^-1))^-1;
    % Z=Z_wf_cable_g;
    Z=Z_wf;
    Z2=Z_wf_cable_g3;
    % Z2=Z_wf;
    clear mag;
    clear pass;
        %MIMO passivity
    pass=mimo_passivity(Z);

    mag=abs(Z2.ResponseData);
    % pass=real(Z(out,in).ResponseData);
    % pass=min(eig_f,[],1);
    im_part=imag(Z(out,in).ResponseData);
    figure(f);

    p2=loglog(fout,squeeze(mag(out,in,:)),'Color',color_plot(i,:));
    % p2=semilogx(fout,squeeze(im_part(out,in,:)),'Color',color_plot(i,:));
    grid on;hold on;
    p2.LineWidth=1.5;
    % [~,neg_ind]=find(pass(1,1,:)<0.00);
    [~,neg_ind]=find(pass<0.00);
    mag_np=NaN(size(squeeze(mag(1,1,:))));
    mag_np(neg_ind)=squeeze(mag(1,1,neg_ind));
    p_np=loglog(fout,mag_np,'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end

    % nexttile(t2,3);
    % p3=semilogx(fout,squeeze(phase(1,1,:)),'Color',color_plot(i,:));
    % p3.LineWidth=1;
    % grid on;hold on;
    % yticks([-180 -90 0 90 180]);
end


% nexttile(t2,1);
line([Z.Frequency(1) Z.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
set(gca,'XTickLabel',[]);
xlim([fout(1), fout(end)]);
% ylim([-200 300]);
ylabel('$\Re\{Z_{pp}(s)\} (\Omega)$','Interpreter','latex');


% nexttile(t2,2);
% ax=gca;
% ax.FontSize=f_size;
% xlim([fout(1), fout(end)]);
% xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');
% 

c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
% c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$P_{ref}$','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);

%% WF impedance with VSM control parameter variation, grid impedance
% Sb_WT=5e6;
% Sb_WF=(12*1)*5e6;
% n_s=1;
% V_wt=900;
% P_pu=1;
% Q_pu=0;
% phi_wt=0;
% D_WW=1;
% % %Rg
% % rg=[0.001, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03,0.04,0.07, 0.1];
% rg=linspace(0.001,0.1,9);
% lg=0.1;
% cvtr.parameters.lg=lg;
% Lg
rg=0.01;
cvtr.parameters.rg=rg;
% lg=[0.001, 0.005, 0.01, 0.015, 0.02,0.05, 0.1,0.2,0.3,0.4];
lg=linspace(0.01, 0.4,12);

p_vary=lg;

f=figure;
% f=gcf;
f.Position=f.Position.*[0.5,0.5,1.5,0.8];
% t2=tiledlayout(2,1);
% t2.TileSpacing='tight';
% nexttile(t2);

out=1;
in=1;
fout=Z_wf.Frequency;

color_plot = colormap(jet(length(p_vary)));
% color_plot = flip(colormap(jet(length(p_vary))));
min_R=0;
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
    [~,neg_ind]=find(pass<0.00);
    mag_np=NaN(size(squeeze(mag(1,1,:))));
    mag_np(neg_ind)=squeeze(mag(1,1,neg_ind));
    % p_np=loglog(fout(neg_ind),squeeze(mag(out,in,neg_ind)),'Color','#515151');
    p_np=loglog(fout,mag_np,'Color','#515151');
    try
        p_np.LineWidth=1.5;
    end
end


f_size=16;

ax=gca;
ax.FontSize=f_size;
xlim([fout(1), fout(end)]);
xlabel('Frequency (Hz)','interpreter','latex');ylabel('$|Z_{pp}(s)| (\Omega)$','Interpreter','latex');


c=colorbar('FontSize',12,'TickLabelInterpreter','latex');
% c.Layout.Tile = 'east';
clim([p_vary(1) p_vary(end)]);
title(c,'$Param$','interpreter','latex');
% exportgraphics(f,[path_for_save,'Z_wf_turbines_vary','.pdf']);

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
% 
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