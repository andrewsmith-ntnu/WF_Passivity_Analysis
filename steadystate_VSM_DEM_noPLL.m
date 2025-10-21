function conv = steadystate_VSM_DEM_noPLL(conv,input)
% Calculate the steady state operating point of a grid-forming VSM DEM converter without PLL

switch nargin
case 1
% conv.ss.idcs0=conv.ss.idcs0;
case 2
% conv.ss.idcs0=input;
end
%% determine the steady state conditions
fsolve_options=optimoptions('fsolve');
fsolve_options.TolFun=1e-10;
fsolve_options.TolX=1e-10;
fsolve_options.Display='off';

sol_ss=fsolve(@(x)steadystate_nonlinear_VSM_DEM_noPLL(x,conv),[1,0]);
conv.ss.DeltaThetavsm0 = sol_ss(1);
conv.ss.iod0 = sol_ss(2);

conv.ss.ioq0 = -1./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(conv.ss.iod0.*conv.ss.Omegag0.*(conv.parameters.lg+conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.Omegag0.^2)+(conv.ss.vnq0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnd0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^2).*cos(conv.ss.DeltaThetavsm0)+( ...
  -conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(conv.ss.DeltaThetavsm0));
conv.ss.voq0 = -1./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(conv.ss.iod0.*conv.ss.Omegag0.*(-conv.parameters.lg.*conv.parameters.rs+conv.parameters.rg.*(conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs)+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.rs.*conv.ss.Omegag0.^2)+conv.parameters.rs.*(conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.ss.vnq0.*(-1+conv.parameters.cf.*conv.parameters.lg.*conv.ss.Omegag0.^2)).*cos( ...
  conv.ss.DeltaThetavsm0)+conv.parameters.rs.*(conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(conv.ss.DeltaThetavsm0));
conv.ss.vod0 = -1./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(conv.ss.iod0.*(-conv.parameters.rg.*(conv.parameters.rg+conv.parameters.rs)-conv.parameters.lg.*(conv.parameters.lg+conv.parameters.ls).*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.^2.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.ls.*conv.ss.Omegag0.^4)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnd0-conv.parameters.lg.*conv.ss.vnq0.*conv.ss.Omegag0+ ...
  conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^3).*cos(conv.ss.DeltaThetavsm0)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnq0+conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0.^2-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^3).*sin( ...
  conv.ss.DeltaThetavsm0));
conv.ss.io0 = conv.ss.iod0+1i.*conv.ss.ioq0;
conv.ss.vo0 = conv.ss.vod0+1i.*conv.ss.voq0;
conv.ss.vir0 = conv.ss.io0.*(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0)+conv.ss.vo0.*(1+conv.parameters.cf.*conv.ss.Omegag0.*(1i.*conv.parameters.rf-conv.parameters.lf.*conv.ss.Omegag0));
conv.ss.il0 = conv.ss.io0+1i.*conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0;
conv.ss.vi0 = conv.ss.io0.*(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0)+conv.ss.vo0.*(1+conv.parameters.cf.*conv.ss.Omegag0.*(1i.*conv.parameters.rf-conv.parameters.lf.*conv.ss.Omegag0));
conv.ss.ilr0 = conv.ss.io0+1i.*conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0;
conv.ss.Gamma0 = 1./conv.parameters.kic.*(conv.ss.io0.*conv.parameters.rf+conv.ss.vo0-conv.parameters.kffv.*conv.ss.vo0+1i.*conv.parameters.cf.*conv.parameters.rf.*conv.ss.vo0.*conv.ss.Omegag0);
conv.ss.Phi0 = conv.ss.vo0;
conv.ss.p0 = real(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.q0 = imag(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.qm0 = imag(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.ve0 = conv.parameters.kdrpq.*conv.ss.qref0+conv.ss.vref0-conv.parameters.kdrpq.*imag(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.Kappa0 = conv.ss.Omegag0;
conv.ss.Omegavsm0 = conv.ss.Omegag0;
conv.ss.vb0 = exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0);

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

% Setting three phase initial conditions
conv.ss.io0a=(cos(conv.ss.DeltaThetavsm0)*real(conv.ss.io0)-sin(conv.ss.DeltaThetavsm0)*imag(conv.ss.io0))*conv.pu.Ib;
conv.ss.io0b=(cos(conv.ss.DeltaThetavsm0-2/3*pi)*real(conv.ss.io0)-sin(conv.ss.DeltaThetavsm0-2/3*pi)*imag(conv.ss.io0))*conv.pu.Ib;
conv.ss.io0c=(cos(conv.ss.DeltaThetavsm0+2/3*pi)*real(conv.ss.io0)-sin(conv.ss.DeltaThetavsm0+2/3*pi)*imag(conv.ss.io0))*conv.pu.Ib;

conv.ss.il0a=(cos(conv.ss.DeltaThetavsm0)*real(conv.ss.il0)-sin(conv.ss.DeltaThetavsm0)*imag(conv.ss.il0))*conv.pu.Ib;
conv.ss.il0b=(cos(conv.ss.DeltaThetavsm0-2/3*pi)*real(conv.ss.il0)-sin(conv.ss.DeltaThetavsm0-2/3*pi)*imag(conv.ss.il0))*conv.pu.Ib;
conv.ss.il0c=(cos(conv.ss.DeltaThetavsm0+2/3*pi)*real(conv.ss.il0)-sin(conv.ss.DeltaThetavsm0+2/3*pi)*imag(conv.ss.il0))*conv.pu.Ib;

conv.ss.vo0a=(cos(conv.ss.DeltaThetavsm0)*real(conv.ss.vo0)-sin(conv.ss.DeltaThetavsm0)*imag(conv.ss.vo0))*conv.pu.Vb;
conv.ss.vo0b=(cos(conv.ss.DeltaThetavsm0-2/3*pi)*real(conv.ss.vo0)-sin(conv.ss.DeltaThetavsm0-2/3*pi)*imag(conv.ss.vo0))*conv.pu.Vb;
conv.ss.vo0c=(cos(conv.ss.DeltaThetavsm0+2/3*pi)*real(conv.ss.vo0)-sin(conv.ss.DeltaThetavsm0+2/3*pi)*imag(conv.ss.vo0))*conv.pu.Vb;