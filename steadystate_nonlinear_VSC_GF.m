function y=steadystate_nonlinear_VSC_GF(x,conv)
% Steady-state nonlinear equations for a grid-following VSC

y= [1./conv.parameters.lg.*conv.pu.Omegab.*(1./x(1).*conv.ss.pref0-x(1).*conv.parameters.rg-x(1).*conv.parameters.lg./conv.ss.pref0.*conv.ss.qref0.*conv.ss.Omegag0-conv.ss.vnd0.*cos(conv.ss.DeltaThetapll0)+conv.ss.vnq0.*sin(conv.ss.DeltaThetapll0)),-1./conv.parameters.lg.*conv.pu.Omegab.*(-x(1)./conv.ss.pref0.*conv.ss.qref0.*conv.parameters.rg+x(1).*conv.parameters.lg.* ...
  conv.ss.Omegag0+conv.ss.vnq0.*cos(conv.ss.DeltaThetapll0)+conv.ss.vnd0.*sin(conv.ss.DeltaThetapll0))];
end