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
% sys = ss(cvtr.A, cvtr.B, cvtr.C, zeros(height(cvtr.C), width(cvtr.B)));
sys = cvtr.sys;

% Select voltage (vo) and current (io) channels in dq ordering
% Note: indices [3 4] and [5 6] depend on how cvtr.C maps states to outputs


inputNames = sys.InputName;
outputNames = sys.OutputName;

% Find correct input
vnd_input = find(contains(inputNames, 'vnd'));
vnq_input = find(contains(inputNames, 'vnq'));

%Find correct outputs
vod_output = find(contains(outputNames, 'vod'));
voq_output = find(contains(outputNames, 'voq'));
iod_output = find(contains(outputNames, 'iod'));
ioq_output = find(contains(outputNames, 'ioq'));
% if any(find(contains(outputNames, 'Thetavsm')))
%     theta_output = find(contains(outputNames, 'Thetavsm')); 
% else
%     theta_output = find(contains(outputNames, 'Theta')); 
% end

% sys_vo_vn = sys([vod_output voq_output theta_output], [vnd_input vnq_input]);
% sys_vo_vn.InputName = {'vnd','vnq'};
% sys_vo_vn.OutputName = {'vod_conv','voq_conv','Theta'};

sys_vo_vn = sys([vod_output voq_output], [vnd_input vnq_input]);
sys_vo_vn.InputName = {'vnd','vnq'};
sys_vo_vn.OutputName = {'vod_conv','voq_conv'};

sys_io_vn = sys([iod_output ioq_output], [vnd_input vnq_input]);
sys_io_vn.InputName = {'vnd','vnq'};
sys_io_vn.OutputName = {'iod','ioq'};

% Build transfer relationships:
sys_vo_io = sys_vo_vn * sys_io_vn^-1;


if isfield(cvtr.ss, 'DeltaThetavsm0')
    theta0 = cvtr.ss.DeltaThetavsm0;
else
    theta0 = cvtr.ss.DeltaThetapll0;
end

% dR = [-sin(theta0), cos(theta0); -cos(theta0), -sin(theta0)];
% 
% D = zeros(2,3);
% D(:,1) = [1;0];  % Coefficient for vnd
% D(:,2) = [0;1];  % Coefficient for vnq
% D(:,3) =-( dR(:,1)*cvtr.ss.vnd0 + dR(:,2)*cvtr.ss.vnq0);  
% 
% 
% sys_conv_vsm=sys_vo_io;
% sys_rot_gframe = ss([], [], [], D);
% sys_rot_gframe.InputName = {'vod_conv', 'voq_conv','Theta'};
% % sys_grid_cframe.InputName = {'vnd', 'vnq'};
% sys_rot_gframe.OutputName = {'vod', 'voq'};
% 
% sys_conv_gframe=connect(sys_conv_vsm, sys_rot_gframe,{'iod', 'ioq'}, {'vod','voq'});
% 
% Z=sys_conv_gframe;
% 


% Impedance Z = -vo/io (following the chosen sign convention)
Z = -sys_vo_io;


% Create FRD object over frequency vector f
Z_frd = frd(Z, f, 'FrequencyUnit', 'Hz');

% Transform each frequency sample from dq to positive/negative sequence
Az = 1/(sqrt(2)) * [1 1i; 1 -1i];
T_theta = [exp(1i*theta0) 0;0 exp(-1i*theta0)];
for i = 1:length(f)
    % apply similarity transform to the 2x2 block
    Z_frd.ResponseData(:,:,i) = Az * Z_frd.ResponseData(:,:,i) * Az^-1;
    % apply angle rotation
    Z_frd.ResponseData(:,:,i) = T_theta * Z_frd.ResponseData(:,:,i) * T_theta^-1;

end

% Rename inputs/outputs for clarity (positive/negative sequence)
Z_frd.InputName = {'i_p','i_n'};
Z_frd.OutputName = {'v_p','v_n'};



