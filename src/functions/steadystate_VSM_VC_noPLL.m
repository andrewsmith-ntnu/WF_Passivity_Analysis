function conv = steadystate_VSM_VC_noPLL(conv,input)
% Calculate the steady state operating point of a grid-forming VSM VC converter without PLL

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

sol_ss=fsolve(@(x)steadystate_nonlinear_VSM_VC_noPLL(x,conv),[0,0],fsolve_options);
conv.ss.ioq0 = sol_ss(1);
conv.ss.DeltaThetavsm0 = sol_ss(2);

conv.ss.iod0 = 1./(conv.parameters.lg+conv.parameters.ls)./conv.ss.Omegag0.*(-conv.ss.ioq0.*(conv.parameters.rg+conv.parameters.rs)-conv.ss.vnq0.*cos(conv.ss.DeltaThetavsm0)+conv.ss.vnd0.*sin(conv.ss.DeltaThetavsm0));
conv.ss.io0 = conv.ss.iod0+1i.*conv.ss.ioq0;
conv.ss.vir0 = (conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0).*(conv.ss.io0-exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0)./(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0).*(-1+conv.parameters.cf.*conv.ss.Omegag0.*((1i*-1).*conv.parameters.rf+conv.parameters.lf.*conv.ss.Omegag0)).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*( ...
  conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)));
conv.ss.vo0 = exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0));
conv.ss.il0 = 1i.*conv.parameters.cf.*exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0).*conv.ss.Omegag0+conv.ss.io0.*(1+conv.parameters.cf.*conv.ss.Omegag0.*(1i.*conv.parameters.rg-conv.parameters.lg.*conv.ss.Omegag0));
conv.ss.vi0 = (conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0).*(conv.ss.io0-exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0)./(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0).*(-1+conv.parameters.cf.*conv.ss.Omegag0.*((1i*-1).*conv.parameters.rf+conv.parameters.lf.*conv.ss.Omegag0)).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*( ...
  conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)));
conv.ss.ilr0 = 1i.*conv.parameters.cf.*exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0).*conv.ss.Omegag0+conv.ss.io0.*(1+conv.parameters.cf.*conv.ss.Omegag0.*(1i.*conv.parameters.rg-conv.parameters.lg.*conv.ss.Omegag0));
conv.ss.Gamma0 = exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0)./conv.parameters.kic.*(-(conv.ss.vnd0+1i.*conv.ss.vnq0).*(-1+conv.parameters.kffv+(1i*-1).*conv.parameters.cf.*conv.parameters.rf.*conv.ss.Omegag0)+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(-(-1+conv.parameters.kffv).*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)+conv.parameters.rf.*(1+conv.parameters.cf.* ...
  conv.ss.Omegag0.*(1i.*conv.parameters.rg-conv.parameters.lg.*conv.ss.Omegag0))));
conv.ss.Phi0 = exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0));
conv.ss.p0 = real(exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)).*conj(conv.ss.io0));
conv.ss.q0 = imag(exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)).*conj(conv.ss.io0));
conv.ss.qm0 = imag(exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)).*conj(conv.ss.io0));
conv.ss.Xi0 = exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0)./conv.parameters.kiv.*(conv.parameters.kpv.*(conv.ss.vnd0+1i.*conv.ss.vnq0)+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*(-conv.parameters.kpv.*(conv.parameters.kdrpq.*conv.ss.qref0+conv.ss.vref0)+conv.ss.io0.*(1-conv.parameters.kffi+conv.parameters.kpv.*(conv.parameters.rg+conv.parameters.rs+1i.*(conv.parameters.lg+conv.parameters.ls).*conv.ss.Omegag0)))+exp(1).^( ...
  1i.*conv.ss.DeltaThetavsm0).*conv.parameters.kdrpq.*conv.parameters.kpv.*imag(exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)).*conj(conv.ss.io0)));
conv.ss.vor0 = conv.parameters.kdrpq.*conv.ss.qref0-conv.ss.io0.*conv.parameters.rs+conv.ss.vref0+(1i*-1).*conv.ss.io0.*conv.parameters.ls.*conv.ss.Omegag0-conv.parameters.kdrpq.*imag(exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)).*conj( ...
  conv.ss.io0));
conv.ss.ve0 = conv.parameters.kdrpq.*conv.ss.qref0+conv.ss.vref0-conv.parameters.kdrpq.*imag(exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0+exp(1).^(1i.*conv.ss.DeltaThetavsm0).*conv.ss.io0.*(conv.parameters.rg+1i.*conv.parameters.lg.*conv.ss.Omegag0)).*conj(conv.ss.io0));
conv.ss.Kappa0 = conv.ss.Omegag0;
conv.ss.Omegavsm0 = conv.ss.Omegag0;
conv.ss.vb0 = exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0);

conv.ss.ild0 = real(conv.ss.il0);
conv.ss.ilq0 = imag(conv.ss.il0);
conv.ss.vid0 = real(conv.ss.vi0);
conv.ss.viq0 = imag(conv.ss.vi0);
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
conv.ss.vird0 = real(conv.ss.vir0);
conv.ss.virq0 = imag(conv.ss.vir0);
conv.ss.vbd0 = real(conv.ss.vb0);
conv.ss.vbq0 = imag(conv.ss.vb0);
conv.ss.vplld0 = real(conv.ss.vpll0);
conv.ss.vpllq0 = imag(conv.ss.vpll0);
conv.ss.Xid0 = real(conv.ss.Xi0);
conv.ss.Xiq0 = imag(conv.ss.Xi0);
conv.ss.ved0 = real(conv.ss.ve0);
conv.ss.veq0 = imag(conv.ss.ve0);
conv.ss.vord0 = real(conv.ss.vor0);
conv.ss.vorq0 = imag(conv.ss.vor0);

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