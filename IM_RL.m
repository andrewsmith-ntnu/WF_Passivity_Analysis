function [Z_frd]=IM_RL(f,f1,R,L)
% IM_RL  Sequence-domain impedance of an RL element
%
% [Z_frd] = IM_RL(f, f1, R, L)
%
% Inputs
%   f  - frequency vector (Hz)
%   f1 - fundamental frequency (Hz)
%   R  - resistance (Ohm)
%   L  - inductance (H)
%
% Outputs
%   Z_frd - FRD object containing 2x2 sequence-domain impedance (pos/neg)

omega = f * 2*pi;
n = 1;
for s = 1i * omega
    % positive and negative sequence impedances (s shifted by +/- f1)
    Zp = R + L * (s + f1*2*pi*1i);
    Zn = R + L * (s - f1*2*pi*1i);
    Z(:,:,n) = [Zp 0; 0 Zn];
    n = n + 1;
end

% return as FRD with Hz units
Z_frd = frd(Z, f, 'FrequencyUnit', 'Hz');