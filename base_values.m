function [cvtr] = base_values(cvtr)
% base_values  Compute per-unit base values for a converter struct
%
% cvtr = base_values(cvtr)
%
% Inputs
%   cvtr - struct with fields 'parameters' containing Vnr, Inr, fnr, etc.
%
% Outputs
%   cvtr - same struct with populated cvtr.pu.* fields used throughout the code

% Base apparent power (kVA)
cvtr.pu.Sb = cvtr.parameters.Vnr * cvtr.parameters.Inr * sqrt(3);
% Base voltage (peak phase voltage)
cvtr.pu.Vb = cvtr.parameters.Vnr * sqrt(2/3);
% Base current (peak)
cvtr.pu.Ib = cvtr.parameters.Inr * sqrt(2);
% Base impedance
cvtr.pu.Zb = cvtr.pu.Vb / cvtr.pu.Ib;
% Base frequency and derived quantities
cvtr.pu.fb = cvtr.parameters.fnr;
cvtr.pu.Omegab = 2*pi*cvtr.pu.fb;
% Base inductance, resistance, capacitance
cvtr.pu.Lb = cvtr.pu.Zb / cvtr.pu.Omegab;
cvtr.pu.Rb = cvtr.pu.Zb;
cvtr.pu.Cb = 1 / (cvtr.pu.Omegab * cvtr.pu.Zb);

% DC-side bases (derived from AC bases)
cvtr.pu.Sbdc = cvtr.pu.Sb;
cvtr.pu.Vbdc = cvtr.pu.Vb * 2;
cvtr.pu.Ibdc = cvtr.pu.Sbdc / cvtr.pu.Vbdc;
cvtr.pu.Rbdc = cvtr.pu.Vbdc / cvtr.pu.Ibdc;
cvtr.pu.Cbdc = 1 / (cvtr.pu.Rbdc * cvtr.pu.Omegab);
cvtr.pu.Lbdc = cvtr.pu.Rbdc / cvtr.pu.Omegab;
end

