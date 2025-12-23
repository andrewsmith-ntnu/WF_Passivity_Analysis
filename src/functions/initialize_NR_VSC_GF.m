function conv = initialize_NR_VSC_GF(conv)
% Initializes the state variables and steady-state values for the Newton-Raphson VSC grid-forming converter model.

conv.ss.ioq0 = -conv.ss.iod0./conv.ss.pref0.*conv.ss.qref0;
conv.ss.voq0 = 0;
conv.ss.vod0 = 1./conv.ss.iod0.*conv.ss.pref0;
conv.ss.io0 = conv.ss.iod0+1i.*conv.ss.ioq0;
conv.ss.vo0 = conv.ss.vod0+1i.*conv.ss.voq0;
conv.ss.il0 = conv.ss.io0+1i.*conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0;
conv.ss.vi0 = conv.ss.vo0+conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0.*(1i.*conv.parameters.rf-conv.parameters.lf.*conv.ss.Omegag0)+conv.ss.io0.*(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0);
conv.ss.ilr0 = conv.ss.io0+1i.*conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0;
conv.ss.Gamma0 = 1./conv.parameters.kic.*(conv.ss.io0.*conv.parameters.rf+conv.ss.vo0-conv.parameters.kffv.*conv.ss.vo0+1i.*conv.parameters.cf.*conv.parameters.rf.*conv.ss.vo0.*conv.ss.Omegag0);
conv.ss.vir0 = conv.ss.vo0+conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0.*(1i.*conv.parameters.rf-conv.parameters.lf.*conv.ss.Omegag0)+conv.ss.io0.*(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0);
conv.ss.Phi0 = conv.ss.vo0;
conv.ss.p0 = real(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.q0 = imag(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.GammaS0 = 1./conv.parameters.kiq.*(conv.ss.io0-conv.parameters.kpp.*conv.ss.pref0+1i.*conv.parameters.kpq.*conv.ss.qref0+1i.*conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0+(1i*-1).*conv.parameters.kpq.*imag(conv.ss.vo0.*conj(conv.ss.io0))+conv.parameters.kpp.*real(conv.ss.vo0.*conj(conv.ss.io0)));
conv.ss.vb0 = exp(1).^(1i.*conv.ss.DeltaThetapll0).*(conv.ss.vnd0+1i.*conv.ss.vnq0);
conv.ss.vpll0 = conv.ss.vo0;
conv.ss.Epsilonpll0 = 1./conv.parameters.kipll.*(conv.ss.Omegag0-conv.parameters.Omegan);
conv.ss.DeltaOmegapll0 = conv.ss.Omegag0-conv.parameters.Omegan;
conv.ss.Omegapll0 = conv.ss.Omegag0;

conv.ss.ild0 = real(conv.ss.il0);
conv.ss.ilq0 = imag(conv.ss.il0);
conv.ss.vid0 = real(conv.ss.vi0);
conv.ss.viq0 = imag(conv.ss.vi0);
conv.ss.vird0 = real(conv.ss.vir0);
conv.ss.virq0 = imag(conv.ss.vir0);
conv.ss.vod0 = real(conv.ss.vo0);
conv.ss.voq0 = imag(conv.ss.vo0);
conv.ss.iod0 = real(conv.ss.io0);
conv.ss.ioq0 = imag(conv.ss.io0);
conv.ss.Gammad0 = real(conv.ss.Gamma0);
conv.ss.Gammaq0 = imag(conv.ss.Gamma0);
conv.ss.Phid0 = real(conv.ss.Phi0);
conv.ss.Phiq0 = imag(conv.ss.Phi0);
conv.ss.ilrd0 = real(conv.ss.ilr0);
conv.ss.ilrq0 = imag(conv.ss.ilr0);
conv.ss.vbd0 = real(conv.ss.vb0);
conv.ss.vbq0 = imag(conv.ss.vb0);
conv.ss.vplld0 = real(conv.ss.vpll0);
conv.ss.vpllq0 = imag(conv.ss.vpll0);
conv.ss.GammaSd0 = real(conv.ss.GammaS0);
conv.ss.GammaSq0 = imag(conv.ss.GammaS0);

conv.x0=[conv.ss.ild0,conv.ss.ilq0,conv.ss.vod0,conv.ss.voq0,conv.ss.iod0,conv.ss.ioq0,conv.ss.Gammad0,conv.ss.Gammaq0,conv.ss.Phid0,conv.ss.Phiq0, ...
  conv.ss.GammaSd0,conv.ss.GammaSq0,conv.ss.vplld0,conv.ss.vpllq0,conv.ss.Epsilonpll0,conv.ss.DeltaThetapll0];
