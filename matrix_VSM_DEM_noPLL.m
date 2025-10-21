function conv = matrix_VSM_DEM_noPLL(conv,conv_str)

%% Calculate the matrixes for the linear model
conv.A=[-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,conv.pu.Omegab.*(conv.ss.Omegag0-conv.ss.Omegavsm0),(-1-conv.parameters.kad+conv.parameters.kffv)./conv.parameters.lf.*conv.pu.Omegab,0,0,0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,0,conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,0,-conv.ss.ilq0.*conv.pu.Omegab,0;
conv.pu.Omegab.*(-conv.ss.Omegag0+conv.ss.Omegavsm0),-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,0,(-1-conv.parameters.kad+conv.parameters.kffv)./conv.parameters.lf.*conv.pu.Omegab,0,0,0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,0,conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,conv.ss.ild0.*conv.pu.Omegab,0;
1./conv.parameters.cf.*conv.pu.Omegab,0,0,conv.pu.Omegab.*conv.ss.Omegag0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0;
0,1./conv.parameters.cf.*conv.pu.Omegab,-conv.pu.Omegab.*conv.ss.Omegag0,0,0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,0;
0,0,1./conv.parameters.lg.*conv.pu.Omegab,0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,conv.pu.Omegab.*conv.ss.Omegag0,0,0,0,0,0,0,0,0,0,1./conv.parameters.lg.*(-conv.ss.vnq0.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0)+conv.ss.vnd0.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0));
0,0,0,1./conv.parameters.lg.*conv.pu.Omegab,-conv.pu.Omegab.*conv.ss.Omegag0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,1./conv.parameters.lg.*conv.pu.Omegab.*(conv.ss.vnd0.*cos(conv.ss.DeltaThetavsm0)+conv.ss.vnq0.*sin(conv.ss.DeltaThetavsm0));
-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
0,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0,0;
0,0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0;
0,0,-conv.ss.ioq0.*conv.parameters.Omegaqf,conv.ss.iod0.*conv.parameters.Omegaqf,conv.ss.voq0.*conv.parameters.Omegaqf,-conv.ss.vod0.*conv.parameters.Omegaqf,0,0,0,0,-conv.parameters.Omegaqf,0,0,0,0,0;
0,0,-1./conv.parameters.ls.*conv.pu.Omegab,0,0,0,0,0,0,0,-conv.parameters.kdrpq./conv.parameters.ls.*conv.pu.Omegab,-1./conv.parameters.ls.*conv.parameters.rs.*conv.pu.Omegab,conv.pu.Omegab.*conv.ss.Omegavsm0,0,conv.ss.ilrq0.*conv.pu.Omegab,0;
0,0,0,-1./conv.parameters.ls.*conv.pu.Omegab,0,0,0,0,0,0,0,-conv.pu.Omegab.*conv.ss.Omegavsm0,-1./conv.parameters.ls.*conv.parameters.rs.*conv.pu.Omegab,0,-conv.ss.ilrd0.*conv.pu.Omegab,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegadf,conv.parameters.Omegadf,0;
0,0,-conv.ss.iod0./conv.parameters.ta,-conv.ss.ioq0./conv.parameters.ta,-1./conv.parameters.ta.*conv.ss.vod0,-1./conv.parameters.ta.*conv.ss.voq0,0,0,0,0,0,0,0,conv.parameters.kd./conv.parameters.ta,-(conv.parameters.kd+conv.parameters.kdrpOmega)./conv.parameters.ta,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,conv.pu.Omegab,0];

conv.B=[0,0,0,0,0,0,conv.ss.ilq0.*conv.pu.Omegab;
0,0,0,0,0,0,-conv.ss.ild0.*conv.pu.Omegab;
0,0,0,0,0,0,conv.ss.voq0.*conv.pu.Omegab;
0,0,0,0,0,0,-conv.ss.vod0.*conv.pu.Omegab;
0,0,-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0),-1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0),0,0,conv.ss.ioq0.*conv.pu.Omegab;
0,0,1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0),-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0),0,0,-conv.ss.iod0.*conv.pu.Omegab;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,conv.parameters.kdrpq./conv.parameters.ls.*conv.pu.Omegab,0,0,1./conv.parameters.ls.*conv.pu.Omegab,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
1./conv.parameters.ta,0,0,0,0,conv.parameters.kdrpOmega./conv.parameters.ta,0;
0,0,0,0,0,0,-conv.pu.Omegab];

conv.C = eye(size(conv.A,1));
conv.sys = ss(conv.A, conv.B, conv.C, 0, 'statename', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_qm') strcat(conv_str,'_ilrd') strcat(conv_str,'_ilrq') strcat(conv_str,'_Kappa') strcat(conv_str,'_Omegavsm') strcat(conv_str,'_DeltaThetavsm') },'inputname', {strcat(conv_str,'_pref') strcat(conv_str,'_qref') strcat(conv_str,'_vnd') strcat(conv_str,'_vnq') strcat(conv_str,'_vref') strcat(conv_str,'_Omegaref') strcat(conv_str,'_Omegag') },'outputname', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_qm') strcat(conv_str,'_ilrd') strcat(conv_str,'_ilrq') strcat(conv_str,'_Kappa') strcat(conv_str,'_Omegavsm') strcat(conv_str,'_DeltaThetavsm') });
end