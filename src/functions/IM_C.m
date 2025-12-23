function [Z_frd]=IM_C(f,f1,C)
% Sequence-domain impedance for a capacitor
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
    % positive and negative sequence impedances
    Zp = 1/(C * (s + f1*2*pi*1i));
    Zn = 1/(C * (s - f1*2*pi*1i));
    Z(:,:,n) = [Zp 0; 0 Zn];
    n = n + 1;
end

% return as FRD with Hz units
Z_frd = frd(Z, f, 'FrequencyUnit', 'Hz');
Z_frd.InputName = {'i_p','i_n'};
Z_frd.OutputName = {'v_p','v_n'};