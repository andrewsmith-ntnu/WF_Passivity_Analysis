function pass = mimo_passivity(Z)
% Simple positive-definiteness check for 2x2 FRD impedance
%
% pass = mimo_passivity(Z)
%
% Inputs
%   Z - 2x2 FRD object representing impedance (Z.ResponseData is 2x2xN)
%
% Outputs
%   pass - vector containing the minimum eigenvalue of the Hermitian part
%          (used as a passivity indicator) for each frequency sample
%
% Notes
%   This routine forms the Hermitian part P = Z + Z^H and computes eigenvalues
%   per frequency. Positive minimum eigenvalue suggests passivity.

% Extract 2x2 FRD blocks
Z_dd = Z.ResponseData(1,1,:);
Z_qd = Z.ResponseData(2,1,:);
Z_dq = Z.ResponseData(1,2,:);
Z_qq = Z.ResponseData(2,2,:);

% Build Hermitian components: a = Z_dd + Z_dd^* etc.
a = Z_dd + conj(Z_dd);
b = Z_qq + conj(Z_qq);
c = Z_dq + conj(Z_qd);

% Construct a PD-like FRD object where each 2x2 sample is [a conj(c); c b]
PD = Z;
PD.ResponseData = [a conj(c); c b];

% For each frequency, compute eigenvalues and collect minimum
for fr = 1:length(PD.Frequency)
    PD_f = PD.ResponseData(:,:,fr);
    sigma = eig(PD_f);
    eig_f(:,fr) = sigma;
end

pass = min(eig_f, [], 1);