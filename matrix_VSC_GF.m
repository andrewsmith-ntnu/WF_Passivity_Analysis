function conv = matrix_VSC_GF(conv,conv_str)

%% Calculate the matrixes for the linear model
conv.A=[-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,0,-(1+conv.parameters.kad-conv.parameters.kffv+conv.ss.iod0.*conv.parameters.kpc.*conv.parameters.kpp)./conv.parameters.lf.*conv.pu.Omegab,-conv.ss.ioq0.*conv.parameters.kpc.*conv.parameters.kpp./conv.parameters.lf.*conv.pu.Omegab,-conv.parameters.kpc.*conv.parameters.kpp./conv.parameters.lf.*conv.ss.vod0.*conv.pu.Omegab,-conv.parameters.kpc.*conv.parameters.kpp./conv.parameters.lf.*conv.ss.voq0.*conv.pu.Omegab,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kiq.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,0,0,0,0;
0,-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,-conv.ss.ioq0.*conv.parameters.kpc.*conv.parameters.kpq./conv.parameters.lf.*conv.pu.Omegab,(-1-conv.parameters.kad+conv.parameters.kffv+conv.ss.iod0.*conv.parameters.kpc.*conv.parameters.kpq)./conv.parameters.lf.*conv.pu.Omegab,conv.parameters.kpc.*conv.parameters.kpq./conv.parameters.lf.*conv.ss.voq0.*conv.pu.Omegab,-conv.parameters.kpc.*conv.parameters.kpq./conv.parameters.lf.*conv.ss.vod0.*conv.pu.Omegab,0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kiq.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,0,0,0;
1./conv.parameters.cf.*conv.pu.Omegab,0,0,conv.pu.Omegab.*(conv.parameters.kipll.*conv.ss.Epsilonpll0+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./conv.ss.vplld0.*conv.ss.vpllq0)),-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,-conv.parameters.kppll.*conv.ss.voq0.*conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,conv.parameters.kppll.*conv.ss.voq0.*conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,conv.parameters.kipll.*conv.ss.voq0.*conv.pu.Omegab,0;
0,1./conv.parameters.cf.*conv.pu.Omegab,-conv.pu.Omegab.*(conv.parameters.kipll.*conv.ss.Epsilonpll0+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./conv.ss.vplld0.*conv.ss.vpllq0)),0,0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,conv.parameters.kppll.*conv.ss.vod0.*conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,-conv.parameters.kppll.*conv.ss.vod0.*conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,-conv.parameters.kipll.*conv.ss.vod0.*conv.pu.Omegab,0;
0,0,1./conv.parameters.lg.*conv.pu.Omegab,0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,conv.pu.Omegab.*(conv.parameters.kipll.*conv.ss.Epsilonpll0+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./conv.ss.vplld0.*conv.ss.vpllq0)),0,0,0,0,0,0,-conv.ss.ioq0.*conv.parameters.kppll.*conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,conv.ss.ioq0.*conv.parameters.kppll.*conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,conv.ss.ioq0.*conv.parameters.kipll.*conv.pu.Omegab,1./conv.parameters.lg.*conv.pu.Omegab.*(conv.ss.vnq0.*cos(conv.ss.DeltaThetapll0)+conv.ss.vnd0.*sin(conv.ss.DeltaThetapll0));
0,0,0,1./conv.parameters.lg.*conv.pu.Omegab,-conv.pu.Omegab.*(conv.parameters.kipll.*conv.ss.Epsilonpll0+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./conv.ss.vplld0.*conv.ss.vpllq0)),-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,0,0,0,0,0,0,conv.ss.iod0.*conv.parameters.kppll.*conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,-conv.ss.iod0.*conv.parameters.kppll.*conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,-conv.ss.iod0.*conv.parameters.kipll.*conv.pu.Omegab,1./conv.parameters.lg.*(-conv.ss.vnd0.*conv.pu.Omegab.*cos(conv.ss.DeltaThetapll0)+conv.ss.vnq0.*conv.pu.Omegab.*sin(conv.ss.DeltaThetapll0));
-1,0,-conv.ss.iod0.*conv.parameters.kpp,-conv.ss.ioq0.*conv.parameters.kpp,-conv.parameters.kpp.*conv.ss.vod0,-conv.parameters.kpp.*conv.ss.voq0,0,0,0,0,conv.parameters.kiq,0,0,0,0,0;
0,-1,-conv.ss.ioq0.*conv.parameters.kpq,conv.ss.iod0.*conv.parameters.kpq,conv.parameters.kpq.*conv.ss.voq0,-conv.parameters.kpq.*conv.ss.vod0,0,0,0,0,0,conv.parameters.kiq,0,0,0,0;
0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0,0;
0,0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0;
0,0,-conv.ss.iod0,-conv.ss.ioq0,-conv.ss.vod0,-conv.ss.voq0,0,0,0,0,0,0,0,0,0,0;
0,0,-conv.ss.ioq0,conv.ss.iod0,conv.ss.voq0,-conv.ss.vod0,0,0,0,0,0,0,0,0,0,0;
0,0,conv.parameters.Omegalppll,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegalppll,0,0,0;
0,0,0,conv.parameters.Omegalppll,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegalppll,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,-conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2),conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2),0,0;
0,0,0,0,0,0,0,0,0,0,0,0,conv.parameters.kppll.*conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,-conv.parameters.kppll.*conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,-conv.parameters.kipll.*conv.pu.Omegab,0];

conv.B=[conv.parameters.kpc.*conv.parameters.kpp./conv.parameters.lf.*conv.pu.Omegab,0,0,0,0;
0,-conv.parameters.kpc.*conv.parameters.kpq./conv.parameters.lf.*conv.pu.Omegab,0,0,0;
0,0,0,0,0;
0,0,0,0,0;
0,0,-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetapll0),1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetapll0),0;
0,0,-1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetapll0),-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetapll0),0;
conv.parameters.kpp,0,0,0,0;
0,-conv.parameters.kpq,0,0,0;
0,0,0,0,0;
0,0,0,0,0;
1,0,0,0,0;
0,-1,0,0,0;
0,0,0,0,0;
0,0,0,0,0;
0,0,0,0,0;
0,0,0,0,conv.pu.Omegab];

conv.C = eye(size(conv.A,1));
conv.sys = ss(conv.A, conv.B, conv.C, 0, 'statename', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_GammaSd') strcat(conv_str,'_GammaSq') strcat(conv_str,'_vplld') strcat(conv_str,'_vpllq') strcat(conv_str,'_Epsilonpll') strcat(conv_str,'_DeltaThetapll') },'inputname', {strcat(conv_str,'_pref') strcat(conv_str,'_qref') strcat(conv_str,'_vnd') strcat(conv_str,'_vnq') strcat(conv_str,'_Omegag') },'outputname', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_GammaSd') strcat(conv_str,'_GammaSq') strcat(conv_str,'_vplld') strcat(conv_str,'_vpllq') strcat(conv_str,'_Epsilonpll') strcat(conv_str,'_DeltaThetapll') });
end