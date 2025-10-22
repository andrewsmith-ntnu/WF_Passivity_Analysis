function [Z_frd]=IM_WT(f,cvtr)
% Compute sequence-domain impedance FRD for a WT converter model
%
% [Z_frd] = IM_WT(f, cvtr)
%
% Inputs
%   f    - frequency vector (Hz)
%   cvtr - converter struct containing state-space matrices A,B,C and parameters
%
% Outputs
%   Z_frd - frequency response data (frd) object of the sequence-domain impedance
%
% Notes
%   This function builds a state-space model from cvtr.A/B/C, selects the
%   relevant voltage/current channels, computes the impedance in dq and
%   transforms it to positive/negative sequence domain.


% Build continuous-time state-space system from provided matrices
sys = ss(cvtr.A, cvtr.B, cvtr.C, zeros(height(cvtr.C), width(cvtr.B)));

% Select voltage (vo) and current (io) channels in dq ordering
% Note: indices [3 4] and [5 6] depend on how cvtr.C maps states to outputs
sys_vo_vn = sys([3 4], [3 4]);
sys_vo_vn.InputName = {'vnd','vnq'};
sys_vo_vn.OutputName = {'vod','voq'};

sys_io_vn = sys([5 6], [3 4]);
sys_io_vn.InputName = {'vnd','vnq'};
sys_io_vn.OutputName = {'iod','ioq'};

% Build transfer relationships: vo/input -> io and vice versa
sys_vo_io = sys_vo_vn * sys_io_vn^-1;
sys_io_vo = sys_io_vn * sys_vo_vn^-1;

% Impedance Z = -vo/io (following the chosen sign convention)
Z = -sys_vo_io;
%Y = -sys_io_vo; % admittance (not used)

% Create FRD object over frequency vector f
Z_frd = frd(Z, f, 'FrequencyUnit', 'Hz');

% Transform each frequency sample from dq to positive/negative sequence
Az = 1/(sqrt(2)) * [1 1i; 1 -1i];
for i = 1:length(f)
    % apply similarity transform to the 2x2 block
    Z_frd.ResponseData(:,:,i) = Az * Z_frd.ResponseData(:,:,i) * Az^-1;
end

% Rename inputs/outputs for clarity (positive/negative sequence)
Z_frd.InputName = {'i_p','i_n'};
Z_frd.OutputName = {'v_p','v_n'};
