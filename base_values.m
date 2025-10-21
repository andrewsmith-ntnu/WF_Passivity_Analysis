function [cvtr] = base_values(cvtr)
% Base values
cvtr.pu.Sb = cvtr.parameters.Vnr*cvtr.parameters.Inr*sqrt(3); % Defining base kVA rating
cvtr.pu.Vb = cvtr.parameters.Vnr*sqrt(2/3); % Base value for voltage defined as peak value of phase voltage
% cvtr.pu.Vb = cvtr.parameters.Vnr*sqrt(1/3); % Base value for voltage defined as peak value of phase voltage
cvtr.pu.Ib = cvtr.parameters.Inr*sqrt(2); % Base value for currents defined as peak value or rated RMS current
% cvtr.pu.Ib = cvtr.parameters.Inr; % Base value for currents defined as peak value or rated RMS current
cvtr.pu.Zb = cvtr.pu.Vb/cvtr.pu.Ib; % Base value for impedance
cvtr.pu.fb = cvtr.parameters.fnr; % Base value for grid frequency
cvtr.pu.Omegab = 2*pi*cvtr.pu.fb;
cvtr.pu.Lb = cvtr.pu.Zb/cvtr.pu.Omegab; % Base for filter inductance
cvtr.pu.Rb = cvtr.pu.Zb;
cvtr.pu.Cb = 1/(cvtr.pu.Omegab*cvtr.pu.Zb);

cvtr.pu.Sbdc = cvtr.pu.Sb;
cvtr.pu.Vbdc = cvtr.pu.Vb*2;
cvtr.pu.Ibdc = cvtr.pu.Sbdc/cvtr.pu.Vbdc;
cvtr.pu.Rbdc = cvtr.pu.Vbdc/cvtr.pu.Ibdc;
cvtr.pu.Cbdc = 1/(cvtr.pu.Rbdc*cvtr.pu.Omegab);
cvtr.pu.Lbdc = cvtr.pu.Rbdc/cvtr.pu.Omegab;
end

