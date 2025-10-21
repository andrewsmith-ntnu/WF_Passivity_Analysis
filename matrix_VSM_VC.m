function conv = matrix_VSM_VC(conv,conv_str)
% Calculate the matrixes for the linear model of the VSM VC converter with PLL.  

%% Calculate the matrixes for the linear model
conv.A=[-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,0,-(1+conv.parameters.kad-conv.parameters.kffv+conv.parameters.kpc.*conv.parameters.kpv)./conv.parameters.lf.*conv.pu.Omegab,-conv.parameters.cf.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab.*conv.ss.Omegag0,conv.parameters.kpc./conv.parameters.lf.*(conv.parameters.kffi-conv.parameters.kpv.*conv.parameters.rs).*conv.pu.Omegab,conv.parameters.kpc.*conv.parameters.kpv./conv.parameters.lf.*conv.parameters.ls.*conv.pu.Omegab.*conv.ss.Omegag0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,-conv.parameters.kdrpq.*conv.parameters.kpc.*conv.parameters.kpv./conv.parameters.lf.*conv.pu.Omegab,conv.parameters.kiv.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,0,0,0,-1./conv.parameters.lf.*(conv.ss.ilq0.*conv.parameters.lf-conv.ss.ioq0.*conv.parameters.kpc.*conv.parameters.kpv.*conv.parameters.ls+conv.parameters.cf.*conv.parameters.kpc.*conv.ss.voq0).*conv.pu.Omegab,0,0;
0,-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,conv.parameters.cf.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab.*conv.ss.Omegag0,-(1+conv.parameters.kad-conv.parameters.kffv+conv.parameters.kpc.*conv.parameters.kpv)./conv.parameters.lf.*conv.pu.Omegab,-conv.parameters.kpc.*conv.parameters.kpv./conv.parameters.lf.*conv.parameters.ls.*conv.pu.Omegab.*conv.ss.Omegag0,conv.parameters.kpc./conv.parameters.lf.*(conv.parameters.kffi-conv.parameters.kpv.*conv.parameters.rs).*conv.pu.Omegab,0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,0,conv.parameters.kiv.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,0,0,1./conv.parameters.lf.*(conv.ss.ild0.*conv.parameters.lf-conv.ss.iod0.*conv.parameters.kpc.*conv.parameters.kpv.*conv.parameters.ls+conv.parameters.cf.*conv.parameters.kpc.*conv.ss.vod0).*conv.pu.Omegab,0,0;
1./conv.parameters.cf.*conv.pu.Omegab,0,0,conv.pu.Omegab.*conv.ss.Omegag0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,1./conv.parameters.cf.*conv.pu.Omegab,-conv.pu.Omegab.*conv.ss.Omegag0,0,0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,1./conv.parameters.lg.*conv.pu.Omegab,0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,conv.pu.Omegab.*conv.ss.Omegag0,0,0,0,0,0,0,0,0,0,0,0,1./conv.parameters.lg.*(-conv.ss.vnq0.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0)+conv.ss.vnd0.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0)),0;
0,0,0,1./conv.parameters.lg.*conv.pu.Omegab,-conv.pu.Omegab.*conv.ss.Omegag0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0,1./conv.parameters.lg.*conv.pu.Omegab.*(conv.ss.vnd0.*cos(conv.ss.DeltaThetavsm0)+conv.ss.vnq0.*sin(conv.ss.DeltaThetavsm0)),0;
-1,0,-conv.parameters.kpv,-conv.parameters.cf.*conv.ss.Omegag0,conv.parameters.kffi-conv.parameters.kpv.*conv.parameters.rs,conv.parameters.kpv.*conv.parameters.ls.*conv.ss.Omegag0,0,0,0,0,-conv.parameters.kdrpq.*conv.parameters.kpv,conv.parameters.kiv,0,0,0,0,conv.ss.ioq0.*conv.parameters.kpv.*conv.parameters.ls-conv.parameters.cf.*conv.ss.voq0,0,0;
0,-1,conv.parameters.cf.*conv.ss.Omegag0,-conv.parameters.kpv,-conv.parameters.kpv.*conv.parameters.ls.*conv.ss.Omegag0,conv.parameters.kffi-conv.parameters.kpv.*conv.parameters.rs,0,0,0,0,0,0,conv.parameters.kiv,0,0,0,-conv.ss.iod0.*conv.parameters.kpv.*conv.parameters.ls+conv.parameters.cf.*conv.ss.vod0,0,0;
0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0,0,0,0,0;
0,0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0,0,0,0;
0,0,-conv.ss.ioq0.*conv.parameters.Omegaqf,conv.ss.iod0.*conv.parameters.Omegaqf,conv.ss.voq0.*conv.parameters.Omegaqf,-conv.ss.vod0.*conv.parameters.Omegaqf,0,0,0,0,-conv.parameters.Omegaqf,0,0,0,0,0,0,0,0;
0,0,-1,0,-conv.parameters.rs,conv.parameters.ls.*conv.ss.Omegag0,0,0,0,0,-conv.parameters.kdrpq,0,0,0,0,0,conv.ss.ioq0.*conv.parameters.ls,0,0;
0,0,0,-1,-conv.parameters.ls.*conv.ss.Omegag0,-conv.parameters.rs,0,0,0,0,0,0,0,0,0,0,-conv.ss.iod0.*conv.parameters.ls,0,0;
0,0,conv.parameters.Omegalppll.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),conv.parameters.Omegalppll.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),0,0,0,0,0,0,0,0,0,-conv.parameters.Omegalppll,0,0,0,-conv.ss.voq0.*conv.parameters.Omegalppll.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)+conv.ss.vod0.*conv.parameters.Omegalppll.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),conv.parameters.Omegalppll.*(conv.ss.voq0.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)-conv.ss.vod0.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0));
0,0,-conv.parameters.Omegalppll.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),conv.parameters.Omegalppll.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),0,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegalppll,0,0,conv.parameters.Omegalppll.*(conv.ss.vod0.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)+conv.ss.voq0.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)),-conv.parameters.Omegalppll.*(conv.ss.vod0.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)+conv.ss.voq0.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0));
0,0,0,0,0,0,0,0,0,0,0,0,0,0,1./conv.ss.vplld0,0,0,0,0;
0,0,-conv.ss.iod0./conv.parameters.ta,-conv.ss.ioq0./conv.parameters.ta,-1./conv.parameters.ta.*conv.ss.vod0,-1./conv.parameters.ta.*conv.ss.voq0,0,0,0,0,0,0,0,0,conv.parameters.kd.*conv.parameters.kppll./conv.parameters.ta./conv.ss.vplld0,conv.parameters.kd.*conv.parameters.kipll./conv.parameters.ta,-(conv.parameters.kd+conv.parameters.kdrpOmega)./conv.parameters.ta,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,conv.pu.Omegab,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,conv.parameters.kppll./conv.ss.vplld0.*conv.pu.Omegab,conv.parameters.kipll.*conv.pu.Omegab,0,0,0];

conv.B=[0,conv.parameters.kdrpq.*conv.parameters.kpc.*conv.parameters.kpv./conv.parameters.lf.*conv.pu.Omegab,0,0,conv.parameters.kpc.*conv.parameters.kpv./conv.parameters.lf.*conv.pu.Omegab,0,conv.ss.ilq0.*conv.pu.Omegab;
0,0,0,0,0,0,-conv.ss.ild0.*conv.pu.Omegab;
0,0,0,0,0,0,conv.ss.voq0.*conv.pu.Omegab;
0,0,0,0,0,0,-conv.ss.vod0.*conv.pu.Omegab;
0,0,-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0),-1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0),0,0,conv.ss.ioq0.*conv.pu.Omegab;
0,0,1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0),-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0),0,0,-conv.ss.iod0.*conv.pu.Omegab;
0,conv.parameters.kdrpq.*conv.parameters.kpv,0,0,conv.parameters.kpv,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,conv.parameters.kdrpq,0,0,1,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
1./conv.parameters.ta,0,0,0,0,conv.parameters.kdrpOmega./conv.parameters.ta,0;
0,0,0,0,0,0,-conv.pu.Omegab;
0,0,0,0,0,0,-conv.pu.Omegab];

conv.C = eye(size(conv.A,1));
conv.sys = ss(conv.A, conv.B, conv.C, 0, 'statename', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_qm') strcat(conv_str,'_Xid') strcat(conv_str,'_Xiq') strcat(conv_str,'_vplld') strcat(conv_str,'_vpllq') strcat(conv_str,'_Epsilonpll') strcat(conv_str,'_Omegavsm') strcat(conv_str,'_DeltaThetavsm') strcat(conv_str,'_DeltaThetapll') },'inputname', {strcat(conv_str,'_pref') strcat(conv_str,'_qref') strcat(conv_str,'_vnd') strcat(conv_str,'_vnq') strcat(conv_str,'_vref') strcat(conv_str,'_Omegaref') strcat(conv_str,'_Omegag') },'outputname', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_qm') strcat(conv_str,'_Xid') strcat(conv_str,'_Xiq') strcat(conv_str,'_vplld') strcat(conv_str,'_vpllq') strcat(conv_str,'_Epsilonpll') strcat(conv_str,'_Omegavsm') strcat(conv_str,'_DeltaThetavsm') strcat(conv_str,'_DeltaThetapll') });
end