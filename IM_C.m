% Im c — auto-generated summary
%
% Auto-generated placeholder. Please replace with a concise description and usage.
function [Z_frd]=IM_C(f,f1,C)
% IM_C  Sequence-domain impedance for a capacitor with frequency shift
%
% [Z_frd] = IM_C(f, f1, C)
%
% Inputs
%   f  - frequency vector (Hz)
%   f1 - fundamental frequency (Hz)
%   C  - capacitance (F)
%
% Outputs
%   Z_frd - frequency response data object (2x2 sequence-domain impedance)

omega = f * 2*pi;
n = 1;
for s = 1i * omega
    % positive and negative sequence impedances (frequency shifted by +/- f1)
    Zp = 1/(C * (s + f1*2*pi*1i));
    Zn = 1/(C * (s - f1*2*pi*1i));
    Z(:,:,n) = [Zp 0; 0 Zn];
    n = n + 1;
end

% return as FRD with Hz units
Z_frd = frd(Z, f, 'FrequencyUnit', 'Hz');