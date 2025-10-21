function [cvtr_ss,Z_cvtr_agg] = impedance_WT_blackbox(freq,f1,V_wt,V_coll,Sb_WT,n_t,ref,rg_eq,lg_eq)
% impedance_WT_blackbox  Build black-box WT aggregated impedance model
%
% [cvtr_ss, Z_cvtr_agg] = impedance_WT_blackbox(freq, f1, V_wt, V_coll, Sb_WT, n_t, ref, rg_eq, lg_eq)
%
% Inputs
%   freq   - frequency vector (Hz)
%   f1     - fundamental frequency (Hz)
%   V_wt   - nominal winding voltage of a WT (V)
%   V_coll - collection voltage used for scaling
%   Sb_WT  - rated apparent power of one WT
%   n_t    - number of turbines in the farm
%   ref    - references struct with fields pref, qref, vref, omegaref
%   rg_eq  - equivalent grid resistance for black-box aggregation
%   lg_eq  - equivalent grid inductance for black-box aggregation
%
% Outputs
%   cvtr_ss     - steady-state converter struct used for black-box model
%   Z_cvtr_agg  - aggregated impedance (sequence-domain FRD)

try
  cvtr.parameters.model_selector = 3;
    
    %% preinitialize variables that may be missing
    cvtr.ss.Kappa0 = 0;
    cvtr.ss.Xi0 = 0;
    cvtr.ss.Epsilonpll0 = 0;
    cvtr.ss.DeltaThetapll0 = 0;
    cvtr.ss.vpll0 = 0;
    cvtr.ss.Xi0 = 0;
    cvtr.ss.vomq0 = 0;
    cvtr.ss.vomd0 = 0;
    cvtr.ss.ilqref0=0;
    cvtr.ss.ildref0=0;
    cvtr.ss.Omegavsm0=0;
    cvtr.ss.GammaSd0=0;
    cvtr.ss.GammaSq0=0;
    
    %% Electrical ratings
    cvtr.parameters.Vnr = V_wt; 
    cvtr.parameters.Vnr_pk = cvtr.parameters.Vnr*sqrt(2)/sqrt(3);
    cvtr.parameters.Inr = Sb_WT/cvtr.parameters.Vnr/sqrt(3);
    cvtr.parameters.fnr = f1;
    cvtr.parameters.Omegan = 1;
    cvtr.parameters.Omegag = 1;
    
    cvtr = base_values(cvtr);
    
    % Filter values (pu)
    cvtr.parameters.lf = .1;
    cvtr.parameters.rf = .005;
    cvtr.parameters.cf = 3e-2;
    
    % Grid values (pu)
    cvtr.parameters.rg = 0.01;
    cvtr.parameters.lg = 0.1;
    
    % Parameters in absolute values
    cvtr.parameters.Lf = cvtr.parameters.lf*cvtr.pu.Lb;
    cvtr.parameters.Rf = cvtr.parameters.rf*cvtr.pu.Zb;
    cvtr.parameters.Cf = cvtr.parameters.cf*cvtr.pu.Cb;
    cvtr.parameters.Lg = cvtr.parameters.lg*cvtr.pu.Lb;
    cvtr.parameters.Rg = cvtr.parameters.rg*cvtr.pu.Zb;
    
    %% Controller parameters
    cvtr.parameters.f_sw = 10e3; % switching frequency
    cvtr.parameters.T_sw = 1/(cvtr.parameters.f_sw);    %switching period
    cvtr.parameters.T_cont = 1/cvtr.parameters.f_sw;     %Control sampling time
    cvtr.parameters.Ts_meas = 1/cvtr.parameters.f_sw;
    cvtr.parameters.Ts = 1e-5;        %Simulation time step
    
    opts.inp.ts=cvtr.parameters.T_cont;
    cvtr.parameters.discrete = 0;     %Discrete or continuous controller
    
    if cvtr.parameters.discrete==0
        cvtr.parameters.rt=cvtr.parameters.Ts;
    else
        cvtr.parameters.rt=cvtr.parameters.T_cont;
    end
    
    
    % Tuning current loop PU
    cvtr.parameters.Tl = cvtr.parameters.lf/(cvtr.parameters.rf*cvtr.pu.Omegab);
    cvtr.parameters.Tic = cvtr.parameters.Tl;
    cvtr.parameters.kpc = cvtr.parameters.Tl*cvtr.parameters.rf/(2*cvtr.parameters.T_sw);
    cvtr.parameters.kic = cvtr.parameters.kpc/cvtr.parameters.Tic;
    
    % Tuning voltage loop PU, VCVSM
    cvtr.parameters.av = 2;
    cvtr.parameters.Tc = cvtr.parameters.cf/cvtr.pu.Omegab;
    cvtr.parameters.Teqc = 2*cvtr.parameters.T_sw;
    cvtr.parameters.Tiv = cvtr.parameters.av^2*cvtr.parameters.Teqc;
    cvtr.parameters.kpv = cvtr.parameters.Tc/(cvtr.parameters.av*cvtr.parameters.Teqc);
    cvtr.parameters.kiv = cvtr.parameters.kpv/cvtr.parameters.Tiv;
    
    % Feedforward terms
    cvtr.parameters.kffv = 1;   %Voltage feedforward
    
    % Tuning pll
    cvtr.parameters.apll = 2;
    cvtr.parameters.Tfpll = 0.002;
    cvtr.parameters.Tipll = cvtr.parameters.apll^2*cvtr.parameters.Tfpll;
    cvtr.parameters.kppll = 1/(cvtr.parameters.apll*2*pi*cvtr.parameters.Tfpll);
    cvtr.parameters.kppll = cvtr.parameters.kppll/cvtr.pu.Omegab;
    cvtr.parameters.kipll = cvtr.parameters.kppll/cvtr.parameters.Tipll;
    cvtr.parameters.Omegalppll = 1/cvtr.parameters.Tfpll;
    
    % Active damping
    cvtr.parameters.Tfad = 0.005;
    cvtr.parameters.Omegaad = 1/cvtr.parameters.Tfad;
    cvtr.parameters.kad = 2.075;
    
    % VSM virtual impedance
    cvtr.parameters.rs = .02;
    cvtr.parameters.ls = .5;
    cvtr.parameters.Omegavo = 500;  %QSEM
    
    % Droop and VSM
    cvtr.parameters.ta = 50;
    cvtr.parameters.kd = 40;
    cvtr.parameters.kdrpOmega = 20;
    cvtr.parameters.kdrpq = 0.05;
    cvtr.parameters.Omegadf = 50;
    cvtr.parameters.Omegaqf = 100;
    
    %Active and reactive power controllers, GFL
    cvtr.parameters.ap = 2;
    cvtr.parameters.Tfp = cvtr.parameters.cf/cvtr.pu.Omegab;
    cvtr.parameters.Teqc = 2*cvtr.parameters.T_sw;
    cvtr.parameters.Tip = cvtr.parameters.ap^2*cvtr.parameters.Teqc;
    cvtr.parameters.kpp = cvtr.parameters.Tfp/(cvtr.parameters.ap*cvtr.parameters.Teqc);
    cvtr.parameters.kip = cvtr.parameters.kpp/cvtr.parameters.Tip;
    cvtr.parameters.kpq=cvtr.parameters.kpp;
    cvtr.parameters.kiq=cvtr.parameters.kip;
    
    
    %% References
    cvtr.ss.pref0 = ref.pref;                  
    cvtr.ss.qref0 = ref.qref;                  
    cvtr.ss.vref0 = ref.vref;                  
    cvtr.ss.Omegaref0 = ref.omegaref;              
    cvtr.ss.vn0 = cvtr.parameters.Vnr_pk/cvtr.pu.Vb;    
    cvtr.ss.vnd0 = cvtr.parameters.Vnr_pk/cvtr.pu.Vb;   
    cvtr.ss.vnq0 = 0;                   
    cvtr.ss.Omegag0 = 1;                        
    cvtr.ss.DeltaThetavsm0 = 0;     
    
    
    cvtr_agg=cvtr;
    
    
    %% Electrical ratings
    Sb_WF=n_t*Sb_WT;
    cvtr_agg.parameters.Inr = Sb_WF/cvtr_agg.parameters.Vnr/sqrt(3);
    
    cvtr_agg = base_values(cvtr_agg);
    
    % 
    cvtr_agg.parameters.rg=rg_eq;
    cvtr_agg.parameters.lg=lg_eq;
    
    % Parameters in absolute values
    cvtr_agg.parameters.Lf = cvtr_agg.parameters.lf*cvtr_agg.pu.Lb;
    cvtr_agg.parameters.Rf = cvtr_agg.parameters.rf*cvtr_agg.pu.Zb;
    cvtr_agg.parameters.Cf = cvtr_agg.parameters.cf*cvtr_agg.pu.Cb;
    cvtr_agg.parameters.Lg = cvtr_agg.parameters.lg*cvtr_agg.pu.Lb;
    cvtr_agg.parameters.Rg = cvtr_agg.parameters.rg*cvtr_agg.pu.Zb;
    
    %% Controller parameters
    
    
    % Tuning current loop PU
    cvtr_agg.parameters.Tl = cvtr_agg.parameters.lf/(cvtr_agg.parameters.rf*cvtr_agg.pu.Omegab);
    cvtr_agg.parameters.Tic = cvtr_agg.parameters.Tl;
    cvtr_agg.parameters.kpc = cvtr_agg.parameters.Tl*cvtr_agg.parameters.rf/(2*cvtr_agg.parameters.T_sw);
    cvtr_agg.parameters.kic = cvtr_agg.parameters.kpc/cvtr_agg.parameters.Tic;
    
    % Tuning voltage loop PU, VCVSM
    cvtr_agg.parameters.av = 2;
    cvtr_agg.parameters.Tc = cvtr_agg.parameters.cf/cvtr_agg.pu.Omegab;
    cvtr_agg.parameters.Teqc = 2*cvtr_agg.parameters.T_sw;
    cvtr_agg.parameters.Tiv = cvtr_agg.parameters.av^2*cvtr_agg.parameters.Teqc;
    cvtr_agg.parameters.kpv = cvtr_agg.parameters.Tc/(cvtr_agg.parameters.av*cvtr_agg.parameters.Teqc);
    cvtr_agg.parameters.kiv = cvtr_agg.parameters.kpv/cvtr_agg.parameters.Tiv;
    
    
    % PQ power loop
    cvtr_agg.parameters.ap = 2;
    cvtr_agg.parameters.Tfp = cvtr_agg.parameters.cf/cvtr_agg.pu.Omegab;
    cvtr_agg.parameters.Teqc = 2*cvtr_agg.parameters.T_sw;
    cvtr_agg.parameters.Tip = cvtr_agg.parameters.ap^2*cvtr_agg.parameters.Teqc;
    cvtr_agg.parameters.kpp = cvtr_agg.parameters.Tfp/(cvtr_agg.parameters.ap*cvtr_agg.parameters.Teqc);
    cvtr_agg.parameters.kip = cvtr_agg.parameters.kpp/cvtr_agg.parameters.Tip;
    cvtr_agg.parameters.kpq = cvtr_agg.parameters.kpp;
    cvtr_agg.parameters.kiq = cvtr_agg.parameters.kip;
                                
    
    
    %%
    switch cvtr_agg.parameters.model_selector
       
        case 3
            
            %% Steady state conditions
            cvtr_agg = steadystate_VSM_QSEM(cvtr_agg);
            
            %% Calculate the matrixes
            cvtr_agg = matrix_VSM_QSEM(cvtr_agg,'VSM');
       
    end
    
    eig_A = eig(cvtr_agg.A);
    
    warning on
    if any(real(eig_A)>0)
        warning('System not stable');
    end
    
    
    
    %% Calculate impedance model
    z_cvtr_agg=IM_WT(freq,cvtr_agg);
    Z_cvtr_agg=z_cvtr_agg*cvtr_agg.pu.Zb*(V_coll/cvtr_agg.parameters.Vnr)^2;
    
    
    %% Diagonalize state-space system
    cvtr_ss=canon(cvtr_agg.sys,'modal');
    cvtr_ss=canon(cvtr_ss,'modal');
catch 
    warning('There was a problem calculating the state-space and impedance model. Check input formats.');
end


%% Functions
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

Z=-sys_vo_io;

Z_frd=frd(Z,f,'FrequencyUnit','Hz');

% Change from dq to Sequence domain 
Az=1/(sqrt(2))*[1 1i;1 -1i];

for i=1:length(f)
    Z_frd.ResponseData(:,:,i)=Az*Z_frd.ResponseData(:,:,i)*Az^-1;
end

Z_frd.InputName={'i_p','i_n'};
Z_frd.OutputName={'v_p','v_n'};
end


end

function conv = steadystate_VSM_QSEM(conv,input)

switch nargin
case 1
% conv.ss.idcs0=conv.ss.idcs0;
case 2
% conv.ss.idcs0=input;
end
%% determine the steady state conditions
fsolve_options=optimoptions('fsolve');
fsolve_options.TolFun=1e-10;
fsolve_options.TolX=1e-10;
fsolve_options.Display='off';

sol_ss=fsolve(@(x)steadystate_nonlinear_VSM_QSEM(x,conv),[1,0,0],fsolve_options);
conv.ss.DeltaThetavsm0 = sol_ss(1);
conv.ss.iod0 = sol_ss(2);
conv.ss.DeltaThetapll0 = sol_ss(3);

conv.ss.ioq0 = -1./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(conv.ss.iod0.*conv.ss.Omegag0.*(conv.parameters.lg+conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.Omegag0.^2)+(conv.ss.vnq0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnd0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^2).*cos(conv.ss.DeltaThetavsm0)+( ...
  -conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(conv.ss.DeltaThetavsm0));
conv.ss.voq0 = -1./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(conv.ss.iod0.*conv.ss.Omegag0.*(-conv.parameters.lg.*conv.parameters.rs+conv.parameters.rg.*(conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs)+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.rs.*conv.ss.Omegag0.^2)+conv.parameters.rs.*(conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.ss.vnq0.*(-1+conv.parameters.cf.*conv.parameters.lg.*conv.ss.Omegag0.^2)).*cos( ...
  conv.ss.DeltaThetavsm0)+conv.parameters.rs.*(conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(conv.ss.DeltaThetavsm0));
conv.ss.vod0 = -1./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(conv.ss.iod0.*(-conv.parameters.rg.*(conv.parameters.rg+conv.parameters.rs)-conv.parameters.lg.*(conv.parameters.lg+conv.parameters.ls).*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.^2.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.ls.*conv.ss.Omegag0.^4)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnd0-conv.parameters.lg.*conv.ss.vnq0.*conv.ss.Omegag0+ ...
  conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^3).*cos(conv.ss.DeltaThetavsm0)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnq0+conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0.^2-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^3).*sin( ...
  conv.ss.DeltaThetavsm0));
conv.ss.io0 = conv.ss.iod0+1i.*conv.ss.ioq0;
conv.ss.vo0 = conv.ss.vod0+1i.*conv.ss.voq0;
conv.ss.vir0 = conv.ss.io0.*(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0)+conv.ss.vo0.*(1+conv.parameters.cf.*conv.ss.Omegag0.*(1i.*conv.parameters.rf-conv.parameters.lf.*conv.ss.Omegag0));
conv.ss.il0 = conv.ss.io0+1i.*conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0;
conv.ss.vi0 = conv.ss.io0.*(conv.parameters.rf+1i.*conv.parameters.lf.*conv.ss.Omegag0)+conv.ss.vo0.*(1+conv.parameters.cf.*conv.ss.Omegag0.*(1i.*conv.parameters.rf-conv.parameters.lf.*conv.ss.Omegag0));
conv.ss.ilr0 = conv.ss.io0+1i.*conv.parameters.cf.*conv.ss.vo0.*conv.ss.Omegag0;
conv.ss.Gamma0 = 1./conv.parameters.kic.*(conv.ss.io0.*conv.parameters.rf+conv.ss.vo0-conv.parameters.kffv.*conv.ss.vo0+1i.*conv.parameters.cf.*conv.parameters.rf.*conv.ss.vo0.*conv.ss.Omegag0);
conv.ss.Phi0 = conv.ss.vo0;
conv.ss.p0 = real(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.q0 = imag(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.qm0 = imag(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.ve0 = conv.parameters.kdrpq.*conv.ss.qref0+conv.ss.vref0-conv.parameters.kdrpq.*imag(conv.ss.vo0.*conj(conv.ss.io0));
conv.ss.Omegavsm0 = conv.ss.Omegag0;
conv.ss.vb0 = exp(1).^((1i*-1).*conv.ss.DeltaThetavsm0).*(conv.ss.vnd0+1i.*conv.ss.vnq0);
conv.ss.vom0 = conv.ss.vo0;
conv.ss.vpll0 = exp(1).^((1i*-1).*(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)).*conv.ss.vo0;
conv.ss.Epsilonpll0 = 1./conv.parameters.kipll.*(conv.ss.Omegag0-conv.parameters.Omegan);
conv.ss.DeltaOmegapll0 = conv.ss.Omegag0-conv.parameters.Omegan;
conv.ss.Omegapll0 = conv.ss.Omegag0;

conv.ss.ild0 = real(conv.ss.il0);
conv.ss.ilq0 = imag(conv.ss.il0);
conv.ss.vid0 = real(conv.ss.vi0);
conv.ss.viq0 = imag(conv.ss.vi0);
conv.ss.vird0 = real(conv.ss.vir0);
conv.ss.virq0 = imag(conv.ss.vir0);
conv.ss.vod0 = real(conv.ss.vo0);
conv.ss.voq0 = imag(conv.ss.vo0);
conv.ss.iod0 = real(conv.ss.io0);
conv.ss.ioq0 = imag(conv.ss.io0);
conv.ss.Gammad0 = real(conv.ss.Gamma0);
conv.ss.Gammaq0 = imag(conv.ss.Gamma0);
conv.ss.Phid0 = real(conv.ss.Phi0);
conv.ss.Phiq0 = imag(conv.ss.Phi0);
conv.ss.ilrd0 = real(conv.ss.ilr0);
conv.ss.ilrq0 = imag(conv.ss.ilr0);
conv.ss.vbd0 = real(conv.ss.vb0);
conv.ss.vbq0 = imag(conv.ss.vb0);
conv.ss.vomd0 = real(conv.ss.vom0);
conv.ss.vomq0 = imag(conv.ss.vom0);
conv.ss.vplld0 = real(conv.ss.vpll0);
conv.ss.vpllq0 = imag(conv.ss.vpll0);

% Setting three phase initial conditions
conv.ss.io0a=(cos(conv.ss.DeltaThetavsm0)*real(conv.ss.io0)-sin(conv.ss.DeltaThetavsm0)*imag(conv.ss.io0))*conv.pu.Ib;
conv.ss.io0b=(cos(conv.ss.DeltaThetavsm0-2/3*pi)*real(conv.ss.io0)-sin(conv.ss.DeltaThetavsm0-2/3*pi)*imag(conv.ss.io0))*conv.pu.Ib;
conv.ss.io0c=(cos(conv.ss.DeltaThetavsm0+2/3*pi)*real(conv.ss.io0)-sin(conv.ss.DeltaThetavsm0+2/3*pi)*imag(conv.ss.io0))*conv.pu.Ib;

conv.ss.il0a=(cos(conv.ss.DeltaThetavsm0)*real(conv.ss.il0)-sin(conv.ss.DeltaThetavsm0)*imag(conv.ss.il0))*conv.pu.Ib;
conv.ss.il0b=(cos(conv.ss.DeltaThetavsm0-2/3*pi)*real(conv.ss.il0)-sin(conv.ss.DeltaThetavsm0-2/3*pi)*imag(conv.ss.il0))*conv.pu.Ib;
conv.ss.il0c=(cos(conv.ss.DeltaThetavsm0+2/3*pi)*real(conv.ss.il0)-sin(conv.ss.DeltaThetavsm0+2/3*pi)*imag(conv.ss.il0))*conv.pu.Ib;

conv.ss.vo0a=(cos(conv.ss.DeltaThetavsm0)*real(conv.ss.vo0)-sin(conv.ss.DeltaThetavsm0)*imag(conv.ss.vo0))*conv.pu.Vb;
conv.ss.vo0b=(cos(conv.ss.DeltaThetavsm0-2/3*pi)*real(conv.ss.vo0)-sin(conv.ss.DeltaThetavsm0-2/3*pi)*imag(conv.ss.vo0))*conv.pu.Vb;
conv.ss.vo0c=(cos(conv.ss.DeltaThetavsm0+2/3*pi)*real(conv.ss.vo0)-sin(conv.ss.DeltaThetavsm0+2/3*pi)*imag(conv.ss.vo0))*conv.pu.Vb;

end


function conv = matrix_VSM_QSEM(conv,conv_str)

%% Calculate the matrixes for the linear model
conv.A=[-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,conv.pu.Omegab.*(conv.ss.Omegag0-conv.ss.Omegavsm0),(-1-conv.parameters.kad+conv.parameters.kffv)./conv.parameters.lf.*conv.pu.Omegab,0,0,0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,0,-conv.parameters.kdrpq.*conv.parameters.kpc.*conv.parameters.rs.*conv.pu.Omegab./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),-1./conv.parameters.lf.*conv.pu.Omegab.*(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2).^(-2).*(conv.ss.ilq0.*conv.parameters.lf.*(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2).^2+conv.parameters.kpc.*conv.parameters.ls.*(conv.parameters.rs.^2.*conv.ss.vomq0+(-2).*conv.parameters.ls.*conv.parameters.rs.*(conv.parameters.kdrpq.*conv.ss.qm0-conv.parameters.kdrpq.*conv.ss.qref0+conv.ss.vomd0-conv.ss.vref0).*conv.ss.Omegavsm0-conv.parameters.ls.^2.*conv.ss.vomq0.*conv.ss.Omegavsm0.^2)),0,-conv.parameters.kpc.*conv.parameters.rs.*conv.pu.Omegab./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),-conv.parameters.kpc.*conv.parameters.ls.*conv.pu.Omegab.*conv.ss.Omegavsm0./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,0,0;
conv.pu.Omegab.*(-conv.ss.Omegag0+conv.ss.Omegavsm0),-1./conv.parameters.lf.*(conv.parameters.kpc+conv.parameters.rf).*conv.pu.Omegab,0,(-1-conv.parameters.kad+conv.parameters.kffv)./conv.parameters.lf.*conv.pu.Omegab,0,0,0,conv.parameters.kic./conv.parameters.lf.*conv.pu.Omegab,0,conv.parameters.kad./conv.parameters.lf.*conv.pu.Omegab,conv.parameters.kdrpq.*conv.parameters.kpc.*conv.parameters.ls.*conv.pu.Omegab.*conv.ss.Omegavsm0./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),1./conv.parameters.lf.*conv.pu.Omegab.*(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2).^(-2).*(conv.parameters.kdrpq.*conv.parameters.kpc.*conv.parameters.ls.*(conv.ss.qm0-conv.ss.qref0).*(conv.parameters.rs.^2-conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2)+conv.ss.ild0.*conv.parameters.lf.*(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2).^2+conv.parameters.kpc.*conv.parameters.ls.*(conv.parameters.rs.^2.*(conv.ss.vomd0-conv.ss.vref0)+2.*conv.parameters.ls.*conv.parameters.rs.*conv.ss.vomq0.*conv.ss.Omegavsm0+conv.parameters.ls.^2.*(-conv.ss.vomd0+conv.ss.vref0).*conv.ss.Omegavsm0.^2)),0,conv.parameters.kpc.*conv.parameters.ls.*conv.pu.Omegab.*conv.ss.Omegavsm0./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),-conv.parameters.kpc.*conv.parameters.rs.*conv.pu.Omegab./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,0,0;
1./conv.parameters.cf.*conv.pu.Omegab,0,0,conv.pu.Omegab.*conv.ss.Omegag0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,1./conv.parameters.cf.*conv.pu.Omegab,-conv.pu.Omegab.*conv.ss.Omegag0,0,0,-1./conv.parameters.cf.*conv.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,1./conv.parameters.lg.*conv.pu.Omegab,0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,conv.pu.Omegab.*conv.ss.Omegag0,0,0,0,0,0,0,1./conv.parameters.lg.*(-conv.ss.vnq0.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0)+conv.ss.vnd0.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0)),0,0,0,0,0,0;
0,0,0,1./conv.parameters.lg.*conv.pu.Omegab,-conv.pu.Omegab.*conv.ss.Omegag0,-1./conv.parameters.lg.*conv.parameters.rg.*conv.pu.Omegab,0,0,0,0,0,0,1./conv.parameters.lg.*conv.pu.Omegab.*(conv.ss.vnd0.*cos(conv.ss.DeltaThetavsm0)+conv.ss.vnq0.*sin(conv.ss.DeltaThetavsm0)),0,0,0,0,0,0;
-1,0,0,0,0,0,0,0,0,0,-conv.parameters.kdrpq.*conv.parameters.rs./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),conv.parameters.ls.*(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2).^(-2).*(-conv.parameters.rs.^2.*conv.ss.vomq0+2.*conv.parameters.ls.*conv.parameters.rs.*(conv.parameters.kdrpq.*(conv.ss.qm0-conv.ss.qref0)+conv.ss.vomd0-conv.ss.vref0).*conv.ss.Omegavsm0+conv.parameters.ls.^2.*conv.ss.vomq0.*conv.ss.Omegavsm0.^2),0,-conv.parameters.rs./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),-conv.parameters.ls.*conv.ss.Omegavsm0./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,0,0;
0,-1,0,0,0,0,0,0,0,0,conv.parameters.kdrpq.*conv.parameters.ls.*conv.ss.Omegavsm0./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),conv.parameters.ls.*(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2).^(-2).*(conv.parameters.rs.^2.*(conv.ss.vomd0-conv.ss.vref0)+2.*conv.parameters.ls.*conv.parameters.rs.*conv.ss.vomq0.*conv.ss.Omegavsm0+conv.parameters.ls.^2.*(-conv.ss.vomd0+conv.ss.vref0).*conv.ss.Omegavsm0.^2+conv.parameters.kdrpq.*(conv.ss.qm0-conv.ss.qref0).*(conv.parameters.rs.^2-conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2)),0,conv.parameters.ls.*conv.ss.Omegavsm0./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),-conv.parameters.rs./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,0,0;
0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0,0,0,0,0;
0,0,0,conv.parameters.Omegaad,0,0,0,0,0,-conv.parameters.Omegaad,0,0,0,0,0,0,0,0,0;
0,0,-conv.ss.ioq0.*conv.parameters.Omegaqf,conv.ss.iod0.*conv.parameters.Omegaqf,conv.ss.voq0.*conv.parameters.Omegaqf,-conv.ss.vod0.*conv.parameters.Omegaqf,0,0,0,0,-conv.parameters.Omegaqf,0,0,0,0,0,0,0,0;
0,0,-conv.ss.iod0./conv.parameters.ta,-conv.ss.ioq0./conv.parameters.ta,-1./conv.parameters.ta.*conv.ss.vod0,-1./conv.parameters.ta.*conv.ss.voq0,0,0,0,0,0,-(conv.parameters.kd+conv.parameters.kdrpOmega)./conv.parameters.ta,0,0,0,-conv.parameters.kd.*conv.parameters.kppll.*conv.ss.vpllq0./(conv.parameters.ta.*conv.ss.vplld0.^2+conv.parameters.ta.*conv.ss.vpllq0.^2),conv.parameters.kd.*conv.parameters.kppll.*conv.ss.vplld0./(conv.parameters.ta.*conv.ss.vplld0.^2+conv.parameters.ta.*conv.ss.vpllq0.^2),conv.parameters.kd.*conv.parameters.kipll./conv.parameters.ta,0;
0,0,0,0,0,0,0,0,0,0,0,conv.pu.Omegab,0,0,0,0,0,0,0;
0,0,conv.parameters.Omegavo,0,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegavo,0,0,0,0,0;
0,0,0,conv.parameters.Omegavo,0,0,0,0,0,0,0,0,0,0,-conv.parameters.Omegavo,0,0,0,0;
0,0,conv.parameters.Omegalppll.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),conv.parameters.Omegalppll.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),0,0,0,0,0,0,0,0,-conv.ss.voq0.*conv.parameters.Omegalppll.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)+conv.ss.vod0.*conv.parameters.Omegalppll.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),0,0,-conv.parameters.Omegalppll,0,0,conv.parameters.Omegalppll.*(conv.ss.voq0.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)-conv.ss.vod0.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0));
0,0,-conv.parameters.Omegalppll.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),conv.parameters.Omegalppll.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0),0,0,0,0,0,0,0,0,conv.parameters.Omegalppll.*(conv.ss.vod0.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)+conv.ss.voq0.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)),0,0,0,-conv.parameters.Omegalppll,0,-conv.parameters.Omegalppll.*(conv.ss.vod0.*cos(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0)+conv.ss.voq0.*sin(conv.ss.DeltaThetapll0-conv.ss.DeltaThetavsm0));
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2),conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2),0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-conv.parameters.kppll.*conv.ss.vpllq0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,conv.parameters.kppll.*conv.ss.vplld0./(conv.ss.vplld0.^2+conv.ss.vpllq0.^2).*conv.pu.Omegab,conv.parameters.kipll.*conv.pu.Omegab,0];

conv.B=[0,conv.parameters.kdrpq.*conv.parameters.kpc.*conv.parameters.rs.*conv.pu.Omegab./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,conv.parameters.kpc.*conv.parameters.rs.*conv.pu.Omegab./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,conv.ss.ilq0.*conv.pu.Omegab;
0,-conv.parameters.kdrpq.*conv.parameters.kpc.*conv.parameters.ls.*conv.pu.Omegab.*conv.ss.Omegavsm0./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,-conv.parameters.kpc.*conv.parameters.ls.*conv.pu.Omegab.*conv.ss.Omegavsm0./(conv.parameters.lf.*conv.parameters.rs.^2+conv.parameters.lf.*conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,-conv.ss.ild0.*conv.pu.Omegab;
0,0,0,0,0,0,conv.ss.voq0.*conv.pu.Omegab;
0,0,0,0,0,0,-conv.ss.vod0.*conv.pu.Omegab;
0,0,-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0),-1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0),0,0,conv.ss.ioq0.*conv.pu.Omegab;
0,0,1./conv.parameters.lg.*conv.pu.Omegab.*sin(conv.ss.DeltaThetavsm0),-1./conv.parameters.lg.*conv.pu.Omegab.*cos(conv.ss.DeltaThetavsm0),0,0,-conv.ss.iod0.*conv.pu.Omegab;
0,conv.parameters.kdrpq.*conv.parameters.rs./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,conv.parameters.rs./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0;
0,-conv.parameters.kdrpq.*conv.parameters.ls.*conv.ss.Omegavsm0./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0,-conv.parameters.ls.*conv.ss.Omegavsm0./(conv.parameters.rs.^2+conv.parameters.ls.^2.*conv.ss.Omegavsm0.^2),0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
1./conv.parameters.ta,0,0,0,0,conv.parameters.kdrpOmega./conv.parameters.ta,0;
0,0,0,0,0,0,-conv.pu.Omegab;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,-conv.pu.Omegab];

conv.C = eye(size(conv.A,1));
conv.sys = ss(conv.A, conv.B, conv.C, 0, 'statename', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_qm') strcat(conv_str,'_Omegavsm') strcat(conv_str,'_DeltaThetavsm') strcat(conv_str,'_vomd') strcat(conv_str,'_vomq') strcat(conv_str,'_vplld') strcat(conv_str,'_vpllq') strcat(conv_str,'_Epsilonpll') strcat(conv_str,'_DeltaThetapll') },'inputname', {strcat(conv_str,'_pref') strcat(conv_str,'_qref') strcat(conv_str,'_vnd') strcat(conv_str,'_vnq') strcat(conv_str,'_vref') strcat(conv_str,'_Omegaref') strcat(conv_str,'_Omegag') },'outputname', {strcat(conv_str,'_ild') strcat(conv_str,'_ilq') strcat(conv_str,'_vod') strcat(conv_str,'_voq') strcat(conv_str,'_iod') strcat(conv_str,'_ioq') strcat(conv_str,'_Gammad') strcat(conv_str,'_Gammaq') strcat(conv_str,'_Phid') strcat(conv_str,'_Phiq') strcat(conv_str,'_qm') strcat(conv_str,'_Omegavsm') strcat(conv_str,'_DeltaThetavsm') strcat(conv_str,'_vomd') strcat(conv_str,'_vomq') strcat(conv_str,'_vplld') strcat(conv_str,'_vpllq') strcat(conv_str,'_Epsilonpll') strcat(conv_str,'_DeltaThetapll') });
end

function y=steadystate_nonlinear_VSM_QSEM(x,conv)
y= [1./conv.parameters.ls.*conv.pu.Omegab.*(-x(2).*conv.parameters.rs+conv.ss.vref0-conv.parameters.cf.*conv.parameters.rs.*conv.ss.Omegag0./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(x(2).*conv.ss.Omegag0.*(-conv.parameters.lg.*conv.parameters.rs+conv.parameters.rg.*(conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs)+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.rs.*conv.ss.Omegag0.^2)+conv.parameters.rs.*(conv.parameters.cf.*conv.parameters.rg.* ...
  conv.ss.vnd0.*conv.ss.Omegag0+conv.ss.vnq0.*(-1+conv.parameters.cf.*conv.parameters.lg.*conv.ss.Omegag0.^2)).*cos(x(1))+conv.parameters.rs.*(conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(x(1)))-conv.parameters.ls.*conv.ss.Omegag0./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+ ...
  conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(x(2).*conv.ss.Omegag0.*(conv.parameters.lg+conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.Omegag0.^2)+(conv.ss.vnq0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnd0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^2).*cos(x(1))+(-conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnq0.* ...
  conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(x(1)))+1./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(x(2).*(-conv.parameters.rg.*(conv.parameters.rg+conv.parameters.rs)-conv.parameters.lg.*(conv.parameters.lg+conv.parameters.ls).*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.^2.*conv.ss.Omegag0.^2+conv.parameters.cf.* ...
  conv.parameters.lg.^2.*conv.parameters.ls.*conv.ss.Omegag0.^4)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnd0-conv.parameters.lg.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^3).*cos(x(1))+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnq0+conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.* ...
  conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0.^2-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^3).*sin(x(1)))-conv.parameters.cf.*conv.parameters.ls.*conv.ss.Omegag0.^2./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(x(2).*(-conv.parameters.rg.*(conv.parameters.rg+conv.parameters.rs)-conv.parameters.lg.*(conv.parameters.lg+conv.parameters.ls).*conv.ss.Omegag0.^2+ ...
  conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.^2.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.ls.*conv.ss.Omegag0.^4)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnd0-conv.parameters.lg.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^3).*cos(x(1))+(-(conv.parameters.rg+conv.parameters.rs).* ...
  conv.ss.vnq0+conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0.^2-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^3).*sin(x(1)))+conv.parameters.kdrpq.*(conv.ss.qref0+x(2)./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(x(2).* ...
  conv.ss.Omegag0.*(-conv.parameters.lg.*conv.parameters.rs+conv.parameters.rg.*(conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs)+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.rs.*conv.ss.Omegag0.^2)+conv.parameters.rs.*(conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.ss.vnq0.*(-1+conv.parameters.cf.*conv.parameters.lg.*conv.ss.Omegag0.^2)).*cos(x(1))+conv.parameters.rs.*(conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0-conv.parameters.cf.* ...
  conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(x(1)))+(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).^(-2).*(x(2).*conv.ss.Omegag0.*(conv.parameters.lg+conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.Omegag0.^2)+(conv.ss.vnq0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnd0.* ...
  conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^2).*cos(x(1))+(-conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(x(1))).*(x(2).*(-conv.parameters.rg.*(conv.parameters.rg+conv.parameters.rs)-conv.parameters.lg.*(conv.parameters.lg+conv.parameters.ls).*conv.ss.Omegag0.^2+ ...
  conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.^2.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.ls.*conv.ss.Omegag0.^4)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnd0-conv.parameters.lg.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^3).*cos(x(1))+(-(conv.parameters.rg+conv.parameters.rs).* ...
  conv.ss.vnq0+conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0.^2-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^3).*sin(x(1))))),1./conv.parameters.ta.*(conv.ss.pref0-conv.parameters.kdrpOmega.*conv.ss.Omegag0+conv.parameters.kdrpOmega.*conv.ss.Omegaref0-(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+ ...
  conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).^(-2).*(x(2).*conv.ss.Omegag0.*(-conv.parameters.lg.*conv.parameters.rs+conv.parameters.rg.*(conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs)+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.rs.*conv.ss.Omegag0.^2)+conv.parameters.rs.*(conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.ss.vnq0.*(-1+conv.parameters.cf.*conv.parameters.lg.*conv.ss.Omegag0.^2)).*cos(x(1))+ ...
  conv.parameters.rs.*(conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(x(1))).*(x(2).*conv.ss.Omegag0.*(conv.parameters.lg+conv.parameters.ls+conv.parameters.cf.*conv.parameters.rg.*conv.parameters.rs-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.Omegag0.^2)+(conv.ss.vnq0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnd0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.ls.* ...
  conv.ss.vnq0.*conv.ss.Omegag0.^2).*cos(x(1))+(-conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rs.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(x(1)))+x(2)./(conv.parameters.rg+conv.parameters.rs-conv.parameters.cf.*(conv.parameters.ls.*conv.parameters.rg+conv.parameters.lg.*conv.parameters.rs).*conv.ss.Omegag0.^2).*(x(2).*(-conv.parameters.rg.*(conv.parameters.rg+ ...
  conv.parameters.rs)-conv.parameters.lg.*(conv.parameters.lg+conv.parameters.ls).*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.^2.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.ls.*conv.ss.Omegag0.^4)+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnd0-conv.parameters.lg.*conv.ss.vnq0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0.^2+conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^3).*cos( ...
  x(1))+(-(conv.parameters.rg+conv.parameters.rs).*conv.ss.vnq0+conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.parameters.cf.*conv.parameters.ls.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0.^2-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^3).*sin(x(1)))),1./(conv.parameters.rs.*(-1+conv.parameters.cf.*conv.parameters.lg.*conv.ss.Omegag0.^2)+conv.parameters.rg.*(-1+conv.parameters.cf.*conv.parameters.ls.* ...
  conv.ss.Omegag0.^2)).*(cos(x(3)-x(1)).*(x(2).*conv.ss.Omegag0.*(conv.parameters.ls.*conv.parameters.rg-conv.parameters.lg.*conv.parameters.rs+conv.parameters.cf.*conv.parameters.rg.^2.*conv.parameters.rs+conv.parameters.cf.*conv.parameters.lg.^2.*conv.parameters.rs.*conv.ss.Omegag0.^2)+conv.parameters.rs.*(conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnd0.*conv.ss.Omegag0+conv.ss.vnq0.*(-1+conv.parameters.cf.*conv.parameters.lg.*conv.ss.Omegag0.^2)) ...
  .*cos(x(1))+conv.parameters.rs.*(conv.ss.vnd0+conv.parameters.cf.*conv.parameters.rg.*conv.ss.vnq0.*conv.ss.Omegag0-conv.parameters.cf.*conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0.^2).*sin(x(1)))+sin(x(3)-x(1)).*(x(2).*(conv.parameters.rg.*conv.parameters.rs+conv.parameters.rg.^2.*(1-conv.parameters.cf.*conv.parameters.ls.* ...
  conv.ss.Omegag0.^2)+conv.parameters.lg.*conv.ss.Omegag0.^2.*(conv.parameters.lg+conv.parameters.ls-conv.parameters.cf.*conv.parameters.lg.*conv.parameters.ls.*conv.ss.Omegag0.^2))+(conv.parameters.rs.*conv.ss.vnd0+conv.parameters.lg.*conv.ss.vnq0.*conv.ss.Omegag0.*(1-conv.parameters.cf.*conv.parameters.ls.*conv.ss.Omegag0.^2)+conv.parameters.rg.*(conv.ss.vnd0-conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnd0.*conv.ss.Omegag0.^2)).*cos(x(1))+(conv.parameters.rs.*conv.ss.vnq0+ ...
  conv.parameters.lg.*conv.ss.vnd0.*conv.ss.Omegag0.*(-1+conv.parameters.cf.*conv.parameters.ls.*conv.ss.Omegag0.^2)+conv.parameters.rg.*(conv.ss.vnq0-conv.parameters.cf.*conv.parameters.ls.*conv.ss.vnq0.*conv.ss.Omegag0.^2)).*sin(x(1))))];
end