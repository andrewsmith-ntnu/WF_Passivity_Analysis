function[y,dy]=function_NR_VSC_GF(x,conv)
% Calculates the Newton-Raphson equations and Jacobian for the VSC grid-forming converter model.

y(1) = -1./conv.parameters.lf.*(x(3).*conv.pu.Omegab+x(3).*conv.parameters.kad.*conv.pu.Omegab-x(9).*conv.parameters.kad.*conv.pu.Omegab-x(3).*conv.parameters.kffv.*conv.pu.Omegab-x(7).*conv.parameters.kic.*conv.pu.Omegab+x(1).*conv.parameters.kpc.*conv.pu.Omegab-x(11).*conv.parameters.kiq.*conv.parameters.kpc.*conv.pu.Omegab+x(3).*x(5).*conv.parameters.kpc.*conv.parameters.kpp.*conv.pu.Omegab+x(4).*x(6).*conv.parameters.kpc.*conv.parameters.kpp.*conv.pu.Omegab-conv.parameters.kpc.* ...
  conv.parameters.kpp.*conv.ss.pref0.*conv.pu.Omegab+x(1).*conv.parameters.rf.*conv.pu.Omegab);
y(2) = -1./conv.parameters.lf.*(x(4).*conv.pu.Omegab-x(10).*conv.parameters.kad.*conv.pu.Omegab+x(4).*conv.parameters.kad.*conv.pu.Omegab-x(4).*conv.parameters.kffv.*conv.pu.Omegab-x(8).*conv.parameters.kic.*conv.pu.Omegab+x(2).*conv.parameters.kpc.*conv.pu.Omegab-x(12).*conv.parameters.kiq.*conv.parameters.kpc.*conv.pu.Omegab-x(4).*x(5).*conv.parameters.kpc.*conv.parameters.kpq.*conv.pu.Omegab+x(3).*x(6).*conv.parameters.kpc.*conv.parameters.kpq.*conv.pu.Omegab+ ...
  conv.parameters.kpc.*conv.parameters.kpq.*conv.ss.qref0.*conv.pu.Omegab+x(2).*conv.parameters.rf.*conv.pu.Omegab);
y(3) = -1./conv.parameters.cf.*(-x(1).*conv.pu.Omegab+x(5).*conv.pu.Omegab-x(15).*x(4).*conv.parameters.cf.*conv.parameters.kipll.*conv.pu.Omegab-x(4).*conv.parameters.cf.*conv.pu.Omegab.*conv.parameters.Omegan-x(4).*conv.parameters.cf.*conv.parameters.kppll.*conv.pu.Omegab.*atan(1./x(13).*x(14)));
y(4) = -1./conv.parameters.cf.*(-x(2).*conv.pu.Omegab+x(6).*conv.pu.Omegab+x(15).*x(3).*conv.parameters.cf.*conv.parameters.kipll.*conv.pu.Omegab+x(3).*conv.parameters.cf.*conv.pu.Omegab.*conv.parameters.Omegan+x(3).*conv.parameters.cf.*conv.parameters.kppll.*conv.pu.Omegab.*atan(1./x(13).*x(14)));
y(5) = -1./conv.parameters.lg.*(-x(3).*conv.pu.Omegab-x(15).*x(6).*conv.parameters.kipll.*conv.parameters.lg.*conv.pu.Omegab+x(5).*conv.parameters.rg.*conv.pu.Omegab-x(6).*conv.parameters.lg.*conv.pu.Omegab.*conv.parameters.Omegan-x(6).*conv.parameters.kppll.*conv.parameters.lg.*conv.pu.Omegab.*atan(1./x(13).*x(14))+conv.ss.vnd0.*conv.pu.Omegab.*cos(x(16))-conv.ss.vnq0.*conv.pu.Omegab.*sin(x(16)));
y(6) = -1./conv.parameters.lg.*(-x(4).*conv.pu.Omegab+x(15).*x(5).*conv.parameters.kipll.*conv.parameters.lg.*conv.pu.Omegab+x(6).*conv.parameters.rg.*conv.pu.Omegab+x(5).*conv.parameters.lg.*conv.pu.Omegab.*conv.parameters.Omegan+x(5).*conv.parameters.kppll.*conv.parameters.lg.*conv.pu.Omegab.*atan(1./x(13).*x(14))+conv.ss.vnq0.*conv.pu.Omegab.*cos(x(16))+conv.ss.vnd0.*conv.pu.Omegab.*sin(x(16)));
y(7) = -x(1)+x(11).*conv.parameters.kiq-x(3).*x(5).*conv.parameters.kpp-x(4).*x(6).*conv.parameters.kpp+conv.parameters.kpp.*conv.ss.pref0;
y(8) = -x(2)+x(12).*conv.parameters.kiq+x(4).*x(5).*conv.parameters.kpq-x(3).*x(6).*conv.parameters.kpq-conv.parameters.kpq.*conv.ss.qref0;
y(9) = x(3).*conv.parameters.Omegaad-x(9).*conv.parameters.Omegaad;
y(10) = -x(10).*conv.parameters.Omegaad+x(4).*conv.parameters.Omegaad;
y(11) = -x(3).*x(5)-x(4).*x(6)+conv.ss.pref0;
y(12) = x(4).*x(5)-x(3).*x(6)-conv.ss.qref0;
y(13) = -x(13).*conv.parameters.Omegalppll+x(3).*conv.parameters.Omegalppll;
y(14) = -x(14).*conv.parameters.Omegalppll+x(4).*conv.parameters.Omegalppll;
y(15) = atan(1./x(13).*x(14));
y(16) = -x(15).*conv.parameters.kipll.*conv.pu.Omegab+conv.pu.Omegab.*conv.ss.Omegag0-conv.pu.Omegab.*conv.parameters.Omegan-conv.parameters.kppll.*conv.pu.Omegab.*atan(1./x(13).*x(14));

y = y';

dy=[-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,0,-(1+conv.parameters.kad-conv.parameters.kffv+x(5).*conv.parameters.kpc.*conv.parameters.kpp)./conv.parameters.lf.*conv.pu.Omegab,-x(6).*conv.parameters.kpc.*conv.parameters.kpp./conv.parameters.lf.*conv.pu.Omegab,-x(3).*conv.parameters.kpc.*conv.parameters.kpp./conv.parameters.lf.*conv.pu.Omegab,-x(4).*conv.parameters.kpc.*conv.parameters.kpp./conv.parameters.lf.*conv.pu.Omegab,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kiq.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,0,0,0,0;
0,-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,-x(6).*conv.parameters.kpc.*conv.parameters.kpq./conv.parameters.lf.*conv.pu.Omegab,(-1-conv.parameters.kad+conv.parameters.kffv+x(5).*conv.parameters.kpc.*conv.parameters.kpq)./conv.parameters.lf.*conv.pu.Omegab,x(4).*conv.parameters.kpc.*conv.parameters.kpq./conv.parameters.lf.*conv.pu.Omegab,-x(3).*conv.parameters.kpc.*conv.parameters.kpq./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kiq.*conv.parameters.kpc./conv.parameters.lf.*conv.pu.Omegab,0,0,0,0;
1./conv.parameters.cf.*conv.pu.Omegab,0,0,conv.pu.Omegab.*(x(15).*conv.parameters.kipll+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./x(13).*x(14))),-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,-x(14)./(x(13).^2+x(14).^2).*x(4).*conv.parameters.kppll.*conv.pu.Omegab,x(13)./(x(13).^2+x(14).^2).*x(4).*conv.parameters.kppll.*conv.pu.Omegab,x(4).*conv.parameters.kipll.*conv.pu.Omegab,0;
0,1./conv.parameters.cf.*conv.pu.Omegab,-conv.pu.Omegab.*(x(15).*conv.parameters.kipll+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./x(13).*x(14))),0,0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,x(14)./(x(13).^2+x(14).^2).*x(3).*conv.parameters.kppll.*conv.pu.Omegab,-x(13)./(x(13).^2+x(14).^2).*x(3).*conv.parameters.kppll.*conv.pu.Omegab,-x(3).*conv.parameters.kipll.*conv.pu.Omegab,0;
0,0,1./conv.parameters.lg.*conv.pu.Omegab,0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,conv.pu.Omegab.*(x(15).*conv.parameters.kipll+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./x(13).*x(14))),0,0,0,0,0,0,-x(14)./(x(13).^2+x(14).^2).*x(6).*conv.parameters.kppll.*conv.pu.Omegab,x(13)./(x(13).^2+x(14).^2).*x(6).*conv.parameters.kppll.*conv.pu.Omegab,x(6).*conv.parameters.kipll.*conv.pu.Omegab,1./conv.parameters.lg.*conv.pu.Omegab.*(conv.ss.vnq0.*cos(x(16))+conv.ss.vnd0.*sin(x(16)));
0,0,0,1./conv.parameters.lg.*conv.pu.Omegab,-conv.pu.Omegab.*(x(15).*conv.parameters.kipll+conv.parameters.Omegan+conv.parameters.kppll.*atan(1./x(13).*x(14))),-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,0,0,0,0,0,0,x(14)./(x(13).^2+x(14).^2).*x(5).*conv.parameters.kppll.*conv.pu.Omegab,-x(13)./(x(13).^2+x(14).^2).*x(5).*conv.parameters.kppll.*conv.pu.Omegab,-x(5).*conv.parameters.kipll.*conv.pu.Omegab,1./conv.parameters.lg.*(-conv.ss.vnd0.*conv.pu.Omegab.*cos(x(16))+conv.ss.vnq0.*conv.pu.Omegab.*sin(x(16)));
-1,0,-x(5).*conv.parameters.kpp,-x(6).*conv.parameters.kpp,-x(3).*conv.parameters.kpp,-x(4).*conv.parameters.kpp,0,0,0,0,conv.parameters.kiq,0,0,0,0,0;
0,-1,-x(6).*conv.parameters.kpq,x(5).*conv.parameters.kpq,x(4).*conv.parameters.kpq,-x(3).*conv.parameters.kpq,0,0,0,0,0,conv.parameters.kiq,0,0,0,0;
0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0,0;
0,0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0;
0,0,-x(5),-x(6),-x(3),-x(4),0,0,0,0,0,0,0,0,0,0;
0,0,-x(6),x(5),x(4),-x(3),0,0,0,0,0,0,0,0,0,0;
0,0,conv.parameters.Omegalppll,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegalppll,0,0,0;
0,0,0,conv.parameters.Omegalppll,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegalppll,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,-x(14)./(x(13).^2+x(14).^2),x(13)./(x(13).^2+x(14).^2),0,0;
0,0,0,0,0,0,0,0,0,0,0,0,x(14)./(x(13).^2+x(14).^2).*conv.parameters.kppll.*conv.pu.Omegab,-x(13)./(x(13).^2+x(14).^2).*conv.parameters.kppll.*conv.pu.Omegab,-conv.parameters.kipll.*conv.pu.Omegab,0];

