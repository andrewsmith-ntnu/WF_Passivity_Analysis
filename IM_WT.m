function [Z_frd]=IM_WT(f,cvtr)
%Calculates the impedance in sequence domain from the state-space models.
%This is specific to the models used, note that the outputs are selected
%manually


sys=ss(cvtr.A,cvtr.B,cvtr.C,zeros(height(cvtr.C),width(cvtr.B)));

sys_vo_vn=sys([3 4],[3 4]);
sys_vo_vn.InputName={'vnd','vnq'};
sys_vo_vn.OutputName={'vod','voq'};

sys_io_vn=sys([5 6],[3 4]);
sys_io_vn.InputName={'vnd','vnq'};
sys_io_vn.OutputName={'iod','ioq'};


sys_vo_io=sys_vo_vn*sys_io_vn^-1;
sys_io_vo=sys_io_vn*sys_vo_vn^-1;


Z=-sys_vo_io;
% Y=-sys_io_vo;

Z_frd=frd(Z,f,'FrequencyUnit','Hz');
% Y_frd=frd(Y,f,'FrequencyUnit','Hz');


%% Change from dq to Sequence domain 
Az=1/(sqrt(2))*[1 1i;1 -1i];

for i=1:length(f)
    Z_frd.ResponseData(:,:,i)=Az*Z_frd.ResponseData(:,:,i)*Az^-1;
end

Z_frd.InputName={'i_p','i_n'};
Z_frd.OutputName={'v_p','v_n'};
