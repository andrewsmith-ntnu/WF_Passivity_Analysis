function [cvtr_ss,Z_cvtr_agg] = impedance_WT_blackbox(freq,f1,V_wt,V_coll,Sb_WT,n_t,ref,rg_eq,lg_eq)
try
    id3244032178.parameters.model_selector=3;
    
    %% preinitialize variables that may be missing
    id3244032178.ss.Kappa0 = 0;
    id3244032178.ss.Xi0 = 0;
    id3244032178.ss.Epsilonpll0 = 0;
    id3244032178.ss.DeltaThetapll0 = 0;
    id3244032178.ss.vpll0 = 0;
    id3244032178.ss.Xi0 = 0;
    id3244032178.ss.vomq0 = 0;
    id3244032178.ss.vomd0 = 0;
    id3244032178.ss.ilqref0=0;
    id3244032178.ss.ildref0=0;
    id3244032178.ss.Omegavsm0=0;
    id3244032178.ss.GammaSd0=0;
    id3244032178.ss.GammaSq0=0;
    
    %% Electrical ratings
    id3244032178.parameters.Vnr = V_wt; 
    id3244032178.parameters.Vnr_pk = id3244032178.parameters.Vnr*sqrt(2)/sqrt(3);
    id3244032178.parameters.Inr = Sb_WT/id3244032178.parameters.Vnr/sqrt(3);
    id3244032178.parameters.fnr = f1;
    id3244032178.parameters.Omegan = 1;
    id3244032178.parameters.Omegag = 1;
    
    id3244032178 = id2901768867(id3244032178);
    
    % Filter values (pu)
    id3244032178.parameters.lf = .1;
    id3244032178.parameters.rf = .005;
    id3244032178.parameters.cf = 3e-2;
    
    % Grid values (pu)
    id3244032178.parameters.rg = 0.01;
    id3244032178.parameters.lg = 0.1;
    
    % Parameters in absolute values
    id3244032178.parameters.Lf = id3244032178.parameters.lf*id3244032178.pu.Lb;
    id3244032178.parameters.Rf = id3244032178.parameters.rf*id3244032178.pu.Zb;
    id3244032178.parameters.Cf = id3244032178.parameters.cf*id3244032178.pu.Cb;
    id3244032178.parameters.Lg = id3244032178.parameters.lg*id3244032178.pu.Lb;
    id3244032178.parameters.Rg = id3244032178.parameters.rg*id3244032178.pu.Zb;
    
    %% Controller parameters
    id3244032178.parameters.f_sw = 10e3; % switching frequency
    id3244032178.parameters.T_sw = 1/(id3244032178.parameters.f_sw);    %switching period
    id3244032178.parameters.T_cont = 1/id3244032178.parameters.f_sw;     %Control sampling time
    id3244032178.parameters.Ts_meas = 1/id3244032178.parameters.f_sw;
    id3244032178.parameters.Ts = 1e-5;        %Simulation time step
    
    id2696738421.inp.ts=id3244032178.parameters.T_cont;
    id3244032178.parameters.discrete = 0;     %Discrete or continuous controller
    
    if id3244032178.parameters.discrete==0
        id3244032178.parameters.rt=id3244032178.parameters.Ts;
    else
        id3244032178.parameters.rt=id3244032178.parameters.T_cont;
    end
    
    
    % Tuning current loop PU
    id3244032178.parameters.Tl = id3244032178.parameters.lf/(id3244032178.parameters.rf*id3244032178.pu.Omegab);
    id3244032178.parameters.Tic = id3244032178.parameters.Tl;
    id3244032178.parameters.kpc = id3244032178.parameters.Tl*id3244032178.parameters.rf/(2*id3244032178.parameters.T_sw);
    id3244032178.parameters.kic = id3244032178.parameters.kpc/id3244032178.parameters.Tic;
    
    % Tuning voltage loop PU, VCVSM
    id3244032178.parameters.av = 2;
    id3244032178.parameters.Tc = id3244032178.parameters.cf/id3244032178.pu.Omegab;
    id3244032178.parameters.Teqc = 2*id3244032178.parameters.T_sw;
    id3244032178.parameters.Tiv = id3244032178.parameters.av^2*id3244032178.parameters.Teqc;
    id3244032178.parameters.kpv = id3244032178.parameters.Tc/(id3244032178.parameters.av*id3244032178.parameters.Teqc);
    id3244032178.parameters.kiv = id3244032178.parameters.kpv/id3244032178.parameters.Tiv;
    
    % Feedforward terms
    id3244032178.parameters.kffv = 1;   %Voltage feedforward
    
    % Tuning pll
    id3244032178.parameters.apll = 2;
    id3244032178.parameters.Tfpll = 0.002;
    id3244032178.parameters.Tipll = id3244032178.parameters.apll^2*id3244032178.parameters.Tfpll;
    id3244032178.parameters.kppll = 1/(id3244032178.parameters.apll*2*pi*id3244032178.parameters.Tfpll);
    id3244032178.parameters.kppll = id3244032178.parameters.kppll/id3244032178.pu.Omegab;
    id3244032178.parameters.kipll = id3244032178.parameters.kppll/id3244032178.parameters.Tipll;
    id3244032178.parameters.Omegalppll = 1/id3244032178.parameters.Tfpll;
    
    % Active damping
    id3244032178.parameters.Tfad = 0.005;
    id3244032178.parameters.Omegaad = 1/id3244032178.parameters.Tfad;
    id3244032178.parameters.kad = 2.075;
    
    % VSM virtual impedance
    id3244032178.parameters.rs = .02;
    id3244032178.parameters.ls = .5;
    id3244032178.parameters.Omegavo = 500;  %QSEM
    
    % Droop and VSM
    id3244032178.parameters.ta = 50;
    id3244032178.parameters.kd = 40;
    id3244032178.parameters.kdrpOmega = 20;
    id3244032178.parameters.kdrpq = 0.05;
    id3244032178.parameters.Omegadf = 50;
    id3244032178.parameters.Omegaqf = 100;
    
    %Active and reactive power controllers, GFL
    id3244032178.parameters.ap = 2;
    id3244032178.parameters.Tfp = id3244032178.parameters.cf/id3244032178.pu.Omegab;
    id3244032178.parameters.Teqc = 2*id3244032178.parameters.T_sw;
    id3244032178.parameters.Tip = id3244032178.parameters.ap^2*id3244032178.parameters.Teqc;
    id3244032178.parameters.kpp = id3244032178.parameters.Tfp/(id3244032178.parameters.ap*id3244032178.parameters.Teqc);
    id3244032178.parameters.kip = id3244032178.parameters.kpp/id3244032178.parameters.Tip;
    id3244032178.parameters.kpq=id3244032178.parameters.kpp;
    id3244032178.parameters.kiq=id3244032178.parameters.kip;
    
    
    %% References
    id3244032178.ss.pref0 = ref.pref;                  
    id3244032178.ss.qref0 = ref.qref;                  
    id3244032178.ss.vref0 = ref.vref;                  
    id3244032178.ss.Omegaref0 = ref.omegaref;              
    id3244032178.ss.vn0 = id3244032178.parameters.Vnr_pk/id3244032178.pu.Vb;    
    id3244032178.ss.vnd0 = id3244032178.parameters.Vnr_pk/id3244032178.pu.Vb;   
    id3244032178.ss.vnq0 = 0;                   
    id3244032178.ss.Omegag0 = 1;                        
    id3244032178.ss.DeltaThetavsm0 = 0;     
    
    
    id2540581075=id3244032178;
    
    
    %% Electrical ratings
    id3409525170=n_t*Sb_WT;
    id2540581075.parameters.Inr = id3409525170/id2540581075.parameters.Vnr/sqrt(3);
    
    id2540581075 = id2901768867(id2540581075);
    
    % 
    id2540581075.parameters.rg=rg_eq;
    id2540581075.parameters.lg=lg_eq;
    
    % Parameters in absolute values
    id2540581075.parameters.Lf = id2540581075.parameters.lf*id2540581075.pu.Lb;
    id2540581075.parameters.Rf = id2540581075.parameters.rf*id2540581075.pu.Zb;
    id2540581075.parameters.Cf = id2540581075.parameters.cf*id2540581075.pu.Cb;
    id2540581075.parameters.Lg = id2540581075.parameters.lg*id2540581075.pu.Lb;
    id2540581075.parameters.Rg = id2540581075.parameters.rg*id2540581075.pu.Zb;
    
    %% Controller parameters
    
    
    % Tuning current loop PU
    id2540581075.parameters.Tl = id2540581075.parameters.lf/(id2540581075.parameters.rf*id2540581075.pu.Omegab);
    id2540581075.parameters.Tic = id2540581075.parameters.Tl;
    id2540581075.parameters.kpc = id2540581075.parameters.Tl*id2540581075.parameters.rf/(2*id2540581075.parameters.T_sw);
    id2540581075.parameters.kic = id2540581075.parameters.kpc/id2540581075.parameters.Tic;
    
    % Tuning voltage loop PU, VCVSM
    id2540581075.parameters.av = 2;
    id2540581075.parameters.Tc = id2540581075.parameters.cf/id2540581075.pu.Omegab;
    id2540581075.parameters.Teqc = 2*id2540581075.parameters.T_sw;
    id2540581075.parameters.Tiv = id2540581075.parameters.av^2*id2540581075.parameters.Teqc;
    id2540581075.parameters.kpv = id2540581075.parameters.Tc/(id2540581075.parameters.av*id2540581075.parameters.Teqc);
    id2540581075.parameters.kiv = id2540581075.parameters.kpv/id2540581075.parameters.Tiv;
    
    
    % PQ power loop
    id2540581075.parameters.ap = 2;
    id2540581075.parameters.Tfp = id2540581075.parameters.cf/id2540581075.pu.Omegab;
    id2540581075.parameters.Teqc = 2*id2540581075.parameters.T_sw;
    id2540581075.parameters.Tip = id2540581075.parameters.ap^2*id2540581075.parameters.Teqc;
    id2540581075.parameters.kpp = id2540581075.parameters.Tfp/(id2540581075.parameters.ap*id2540581075.parameters.Teqc);
    id2540581075.parameters.kip = id2540581075.parameters.kpp/id2540581075.parameters.Tip;
    id2540581075.parameters.kpq = id2540581075.parameters.kpp;
    id2540581075.parameters.kiq = id2540581075.parameters.kip;
                                
    
    
    %%
    switch id2540581075.parameters.model_selector
       
        case 3
            
            %% Steady state conditions
            id2540581075 = id1108169701(id2540581075);
            
            %% Calculate the matrixes
            id2540581075 = id3381706392(id2540581075,'VSM');
       
    end
    
    id162062399 = eig(id2540581075.A);
    
    warning on
    if any(real(id162062399)>0)
        warning('System not stable');
    end
    
    
    
    %% Calculate impedance model
    id1280780142=id3638149721(freq,id2540581075);
    Z_cvtr_agg=id1280780142*id2540581075.pu.Zb*(V_coll/id2540581075.parameters.Vnr)^2;
    
    
    %% Diagonalize state-space system
    cvtr_ss=canon(id2540581075.sys,'modal');
    cvtr_ss=canon(cvtr_ss,'modal');
catch 
    warning('There was a problem calculating the state-space and impedance model. Check input formats.');
end


%% Functions
function [id1293877760] = id2901768867(id1293877760)
% Base values
id1293877760.pu.Sb = id1293877760.parameters.Vnr*id1293877760.parameters.Inr*sqrt(3); % Defining base kVA rating
id1293877760.pu.Vb = id1293877760.parameters.Vnr*sqrt(2/3); % Base value for voltage defined as peak value of phase voltage
% cvtr.pu.Vb = cvtr.parameters.Vnr*sqrt(1/3); % Base value for voltage defined as peak value of phase voltage
id1293877760.pu.Ib = id1293877760.parameters.Inr*sqrt(2); % Base value for currents defined as peak value or rated RMS current
% cvtr.pu.Ib = cvtr.parameters.Inr; % Base value for currents defined as peak value or rated RMS current
id1293877760.pu.Zb = id1293877760.pu.Vb/id1293877760.pu.Ib; % Base value for impedance
id1293877760.pu.fb = id1293877760.parameters.fnr; % Base value for grid frequency
id1293877760.pu.Omegab = 2*pi*id1293877760.pu.fb;
id1293877760.pu.Lb = id1293877760.pu.Zb/id1293877760.pu.Omegab; % Base for filter inductance
id1293877760.pu.Rb = id1293877760.pu.Zb;
id1293877760.pu.Cb = 1/(id1293877760.pu.Omegab*id1293877760.pu.Zb);

id1293877760.pu.Sbdc = id1293877760.pu.Sb;
id1293877760.pu.Vbdc = id1293877760.pu.Vb*2;
id1293877760.pu.Ibdc = id1293877760.pu.Sbdc/id1293877760.pu.Vbdc;
id1293877760.pu.Rbdc = id1293877760.pu.Vbdc/id1293877760.pu.Ibdc;
id1293877760.pu.Cbdc = 1/(id1293877760.pu.Rbdc*id1293877760.pu.Omegab);
id1293877760.pu.Lbdc = id1293877760.pu.Rbdc/id1293877760.pu.Omegab;
end


function [id1680148142]=id3638149721(id3399733779,id2522220737)
%Calculates the impedance in sequence domain from the state-space models.
%This is specific to the models used, note that the outputs are selected
%manually

id2331200361=ss(id2522220737.A,id2522220737.B,id2522220737.C,zeros(height(id2522220737.C),width(id2522220737.B)));

id1651967783=id2331200361([3 4],[3 4]);
id1651967783.InputName={'vnd','vnq'};
id1651967783.OutputName={'vod','voq'};

id2921012062=id2331200361([5 6],[3 4]);
id2921012062.InputName={'vnd','vnq'};
id2921012062.OutputName={'iod','ioq'};

id4247374578=id1651967783*id2921012062^-1;

id2231948278=-id4247374578;

id1680148142=frd(id2231948278,id3399733779,'FrequencyUnit','Hz');

% Change from dq to Sequence domain 
id2944043053=1/(sqrt(2))*[1 1i;1 -1i];

for id1456129730=1:length(id3399733779)
    id1680148142.ResponseData(:,:,id1456129730)=id2944043053*id1680148142.ResponseData(:,:,id1456129730)*id2944043053^-1;
end

id1680148142.InputName={'i_p','i_n'};
id1680148142.OutputName={'v_p','v_n'};
end


end

function id1010998369 = id1108169701(id1010998369,id246673237)

switch nargin
case 1
% conv.ss.idcs0=conv.ss.idcs0;
case 2
% conv.ss.idcs0=input;
end
%% determine the steady state conditions
id1070117844=optimoptions('fsolve');
id1070117844.TolFun=1e-10;
id1070117844.TolX=1e-10;
id1070117844.Display='off';

id970760690=fsolve(@(id4084569904)id3325638676(id4084569904,id1010998369),[1,0,0],id1070117844);
id1010998369.ss.DeltaThetavsm0 = id970760690(1);
id1010998369.ss.iod0 = id970760690(2);
id1010998369.ss.DeltaThetapll0 = id970760690(3);

id1010998369.ss.ioq0 = -1./(id1010998369.parameters.rg+id1010998369.parameters.rs-id1010998369.parameters.cf.*(id1010998369.parameters.ls.*id1010998369.parameters.rg+id1010998369.parameters.lg.*id1010998369.parameters.rs).*id1010998369.ss.Omegag0.^2).*(id1010998369.ss.iod0.*id1010998369.ss.Omegag0.*(id1010998369.parameters.lg+id1010998369.parameters.ls+id1010998369.parameters.cf.*id1010998369.parameters.rg.*id1010998369.parameters.rs-id1010998369.parameters.cf.*id1010998369.parameters.lg.*id1010998369.parameters.ls.*id1010998369.ss.Omegag0.^2)+(id1010998369.ss.vnq0+id1010998369.parameters.cf.*id1010998369.parameters.rs.*id1010998369.ss.vnd0.*id1010998369.ss.Omegag0-id1010998369.parameters.cf.*id1010998369.parameters.ls.*id1010998369.ss.vnq0.*id1010998369.ss.Omegag0.^2).*cos(id1010998369.ss.DeltaThetavsm0)+( ...
  -id1010998369.ss.vnd0+id1010998369.parameters.cf.*id1010998369.parameters.rs.*id1010998369.ss.vnq0.*id1010998369.ss.Omegag0+id1010998369.parameters.cf.*id1010998369.parameters.ls.*id1010998369.ss.vnd0.*id1010998369.ss.Omegag0.^2).*sin(id1010998369.ss.DeltaThetavsm0));
id1010998369.ss.voq0 = -1./(id1010998369.parameters.rg+id1010998369.parameters.rs-id1010998369.parameters.cf.*(id1010998369.parameters.ls.*id1010998369.parameters.rg+id1010998369.parameters.lg.*id1010998369.parameters.rs).*id1010998369.ss.Omegag0.^2).*(id1010998369.ss.iod0.*id1010998369.ss.Omegag0.*(-id1010998369.parameters.lg.*id1010998369.parameters.rs+id1010998369.parameters.rg.*(id1010998369.parameters.ls+id1010998369.parameters.cf.*id1010998369.parameters.rg.*id1010998369.parameters.rs)+id1010998369.parameters.cf.*id1010998369.parameters.lg.^2.*id1010998369.parameters.rs.*id1010998369.ss.Omegag0.^2)+id1010998369.parameters.rs.*(id1010998369.parameters.cf.*id1010998369.parameters.rg.*id1010998369.ss.vnd0.*id1010998369.ss.Omegag0+id1010998369.ss.vnq0.*(-1+id1010998369.parameters.cf.*id1010998369.parameters.lg.*id1010998369.ss.Omegag0.^2)).*cos( ...
  id1010998369.ss.DeltaThetavsm0)+id1010998369.parameters.rs.*(id1010998369.ss.vnd0+id1010998369.parameters.cf.*id1010998369.parameters.rg.*id1010998369.ss.vnq0.*id1010998369.ss.Omegag0-id1010998369.parameters.cf.*id1010998369.parameters.lg.*id1010998369.ss.vnd0.*id1010998369.ss.Omegag0.^2).*sin(id1010998369.ss.DeltaThetavsm0));
id1010998369.ss.vod0 = -1./(id1010998369.parameters.rg+id1010998369.parameters.rs-id1010998369.parameters.cf.*(id1010998369.parameters.ls.*id1010998369.parameters.rg+id1010998369.parameters.lg.*id1010998369.parameters.rs).*id1010998369.ss.Omegag0.^2).*(id1010998369.ss.iod0.*(-id1010998369.parameters.rg.*(id1010998369.parameters.rg+id1010998369.parameters.rs)-id1010998369.parameters.lg.*(id1010998369.parameters.lg+id1010998369.parameters.ls).*id1010998369.ss.Omegag0.^2+id1010998369.parameters.cf.*id1010998369.parameters.ls.*id1010998369.parameters.rg.^2.*id1010998369.ss.Omegag0.^2+id1010998369.parameters.cf.*id1010998369.parameters.lg.^2.*id1010998369.parameters.ls.*id1010998369.ss.Omegag0.^4)+(-(id1010998369.parameters.rg+id1010998369.parameters.rs).*id1010998369.ss.vnd0-id1010998369.parameters.lg.*id1010998369.ss.vnq0.*id1010998369.ss.Omegag0+ ...
  id1010998369.parameters.cf.*id1010998369.parameters.ls.*id1010998369.parameters.rg.*id1010998369.ss.vnd0.*id1010998369.ss.Omegag0.^2+id1010998369.parameters.cf.*id1010998369.parameters.lg.*id1010998369.parameters.ls.*id1010998369.ss.vnq0.*id1010998369.ss.Omegag0.^3).*cos(id1010998369.ss.DeltaThetavsm0)+(-(id1010998369.parameters.rg+id1010998369.parameters.rs).*id1010998369.ss.vnq0+id1010998369.parameters.lg.*id1010998369.ss.vnd0.*id1010998369.ss.Omegag0+id1010998369.parameters.cf.*id1010998369.parameters.ls.*id1010998369.parameters.rg.*id1010998369.ss.vnq0.*id1010998369.ss.Omegag0.^2-id1010998369.parameters.cf.*id1010998369.parameters.lg.*id1010998369.parameters.ls.*id1010998369.ss.vnd0.*id1010998369.ss.Omegag0.^3).*sin( ...
  id1010998369.ss.DeltaThetavsm0));
id1010998369.ss.io0 = id1010998369.ss.iod0+1i.*id1010998369.ss.ioq0;
id1010998369.ss.vo0 = id1010998369.ss.vod0+1i.*id1010998369.ss.voq0;
id1010998369.ss.vir0 = id1010998369.ss.io0.*(id1010998369.parameters.rf+1i.*id1010998369.parameters.lf.*id1010998369.ss.Omegag0)+id1010998369.ss.vo0.*(1+id1010998369.parameters.cf.*id1010998369.ss.Omegag0.*(1i.*id1010998369.parameters.rf-id1010998369.parameters.lf.*id1010998369.ss.Omegag0));
id1010998369.ss.il0 = id1010998369.ss.io0+1i.*id1010998369.parameters.cf.*id1010998369.ss.vo0.*id1010998369.ss.Omegag0;
id1010998369.ss.vi0 = id1010998369.ss.io0.*(id1010998369.parameters.rf+1i.*id1010998369.parameters.lf.*id1010998369.ss.Omegag0)+id1010998369.ss.vo0.*(1+id1010998369.parameters.cf.*id1010998369.ss.Omegag0.*(1i.*id1010998369.parameters.rf-id1010998369.parameters.lf.*id1010998369.ss.Omegag0));
id1010998369.ss.ilr0 = id1010998369.ss.io0+1i.*id1010998369.parameters.cf.*id1010998369.ss.vo0.*id1010998369.ss.Omegag0;
id1010998369.ss.Gamma0 = 1./id1010998369.parameters.kic.*(id1010998369.ss.io0.*id1010998369.parameters.rf+id1010998369.ss.vo0-id1010998369.parameters.kffv.*id1010998369.ss.vo0+1i.*id1010998369.parameters.cf.*id1010998369.parameters.rf.*id1010998369.ss.vo0.*id1010998369.ss.Omegag0);
id1010998369.ss.Phi0 = id1010998369.ss.vo0;
id1010998369.ss.p0 = real(id1010998369.ss.vo0.*conj(id1010998369.ss.io0));
id1010998369.ss.q0 = imag(id1010998369.ss.vo0.*conj(id1010998369.ss.io0));
id1010998369.ss.qm0 = imag(id1010998369.ss.vo0.*conj(id1010998369.ss.io0));
id1010998369.ss.ve0 = id1010998369.parameters.kdrpq.*id1010998369.ss.qref0+id1010998369.ss.vref0-id1010998369.parameters.kdrpq.*imag(id1010998369.ss.vo0.*conj(id1010998369.ss.io0));
id1010998369.ss.Omegavsm0 = id1010998369.ss.Omegag0;
id1010998369.ss.vb0 = exp(1).^((1i*-1).*id1010998369.ss.DeltaThetavsm0).*(id1010998369.ss.vnd0+1i.*id1010998369.ss.vnq0);
id1010998369.ss.vom0 = id1010998369.ss.vo0;
id1010998369.ss.vpll0 = exp(1).^((1i*-1).*(id1010998369.ss.DeltaThetapll0-id1010998369.ss.DeltaThetavsm0)).*id1010998369.ss.vo0;
id1010998369.ss.Epsilonpll0 = 1./id1010998369.parameters.kipll.*(id1010998369.ss.Omegag0-id1010998369.parameters.Omegan);
id1010998369.ss.DeltaOmegapll0 = id1010998369.ss.Omegag0-id1010998369.parameters.Omegan;
id1010998369.ss.Omegapll0 = id1010998369.ss.Omegag0;

id1010998369.ss.ild0 = real(id1010998369.ss.il0);
id1010998369.ss.ilq0 = imag(id1010998369.ss.il0);
id1010998369.ss.vid0 = real(id1010998369.ss.vi0);
id1010998369.ss.viq0 = imag(id1010998369.ss.vi0);
id1010998369.ss.vird0 = real(id1010998369.ss.vir0);
id1010998369.ss.virq0 = imag(id1010998369.ss.vir0);
id1010998369.ss.vod0 = real(id1010998369.ss.vo0);
id1010998369.ss.voq0 = imag(id1010998369.ss.vo0);
id1010998369.ss.iod0 = real(id1010998369.ss.io0);
id1010998369.ss.ioq0 = imag(id1010998369.ss.io0);
id1010998369.ss.Gammad0 = real(id1010998369.ss.Gamma0);
id1010998369.ss.Gammaq0 = imag(id1010998369.ss.Gamma0);
id1010998369.ss.Phid0 = real(id1010998369.ss.Phi0);
id1010998369.ss.Phiq0 = imag(id1010998369.ss.Phi0);
id1010998369.ss.ilrd0 = real(id1010998369.ss.ilr0);
id1010998369.ss.ilrq0 = imag(id1010998369.ss.ilr0);
id1010998369.ss.vbd0 = real(id1010998369.ss.vb0);
id1010998369.ss.vbq0 = imag(id1010998369.ss.vb0);
id1010998369.ss.vomd0 = real(id1010998369.ss.vom0);
id1010998369.ss.vomq0 = imag(id1010998369.ss.vom0);
id1010998369.ss.vplld0 = real(id1010998369.ss.vpll0);
id1010998369.ss.vpllq0 = imag(id1010998369.ss.vpll0);

% Setting three phase initial conditions
id1010998369.ss.io0a=(cos(id1010998369.ss.DeltaThetavsm0)*real(id1010998369.ss.io0)-sin(id1010998369.ss.DeltaThetavsm0)*imag(id1010998369.ss.io0))*id1010998369.pu.Ib;
id1010998369.ss.io0b=(cos(id1010998369.ss.DeltaThetavsm0-2/3*pi)*real(id1010998369.ss.io0)-sin(id1010998369.ss.DeltaThetavsm0-2/3*pi)*imag(id1010998369.ss.io0))*id1010998369.pu.Ib;
id1010998369.ss.io0c=(cos(id1010998369.ss.DeltaThetavsm0+2/3*pi)*real(id1010998369.ss.io0)-sin(id1010998369.ss.DeltaThetavsm0+2/3*pi)*imag(id1010998369.ss.io0))*id1010998369.pu.Ib;

id1010998369.ss.il0a=(cos(id1010998369.ss.DeltaThetavsm0)*real(id1010998369.ss.il0)-sin(id1010998369.ss.DeltaThetavsm0)*imag(id1010998369.ss.il0))*id1010998369.pu.Ib;
id1010998369.ss.il0b=(cos(id1010998369.ss.DeltaThetavsm0-2/3*pi)*real(id1010998369.ss.il0)-sin(id1010998369.ss.DeltaThetavsm0-2/3*pi)*imag(id1010998369.ss.il0))*id1010998369.pu.Ib;
id1010998369.ss.il0c=(cos(id1010998369.ss.DeltaThetavsm0+2/3*pi)*real(id1010998369.ss.il0)-sin(id1010998369.ss.DeltaThetavsm0+2/3*pi)*imag(id1010998369.ss.il0))*id1010998369.pu.Ib;

id1010998369.ss.vo0a=(cos(id1010998369.ss.DeltaThetavsm0)*real(id1010998369.ss.vo0)-sin(id1010998369.ss.DeltaThetavsm0)*imag(id1010998369.ss.vo0))*id1010998369.pu.Vb;
id1010998369.ss.vo0b=(cos(id1010998369.ss.DeltaThetavsm0-2/3*pi)*real(id1010998369.ss.vo0)-sin(id1010998369.ss.DeltaThetavsm0-2/3*pi)*imag(id1010998369.ss.vo0))*id1010998369.pu.Vb;
id1010998369.ss.vo0c=(cos(id1010998369.ss.DeltaThetavsm0+2/3*pi)*real(id1010998369.ss.vo0)-sin(id1010998369.ss.DeltaThetavsm0+2/3*pi)*imag(id1010998369.ss.vo0))*id1010998369.pu.Vb;

end


function id2489385430 = id3381706392(id2489385430,id1947728959)

%% Calculate the matrixes for the linear model
id2489385430.A=[-1./id2489385430.parameters.lf.*(id2489385430.parameters.kpc+id2489385430.parameters.rf).*id2489385430.pu.Omegab,id2489385430.pu.Omegab.*(id2489385430.ss.Omegag0-id2489385430.ss.Omegavsm0),(-1-id2489385430.parameters.kad+id2489385430.parameters.kffv)./id2489385430.parameters.lf.*id2489385430.pu.Omegab,0,0,0,id2489385430.parameters.kic./id2489385430.parameters.lf.*id2489385430.pu.Omegab,0,id2489385430.parameters.kad./id2489385430.parameters.lf.*id2489385430.pu.Omegab,0,-id2489385430.parameters.kdrpq.*id2489385430.parameters.kpc.*id2489385430.parameters.rs.*id2489385430.pu.Omegab./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),-1./id2489385430.parameters.lf.*id2489385430.pu.Omegab.*(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2).^(-2).*(id2489385430.ss.ilq0.*id2489385430.parameters.lf.*(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2).^2+id2489385430.parameters.kpc.*id2489385430.parameters.ls.*(id2489385430.parameters.rs.^2.*id2489385430.ss.vomq0+(-2).*id2489385430.parameters.ls.*id2489385430.parameters.rs.*(id2489385430.parameters.kdrpq.*id2489385430.ss.qm0-id2489385430.parameters.kdrpq.*id2489385430.ss.qref0+id2489385430.ss.vomd0-id2489385430.ss.vref0).*id2489385430.ss.Omegavsm0-id2489385430.parameters.ls.^2.*id2489385430.ss.vomq0.*id2489385430.ss.Omegavsm0.^2)),0,-id2489385430.parameters.kpc.*id2489385430.parameters.rs.*id2489385430.pu.Omegab./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),-id2489385430.parameters.kpc.*id2489385430.parameters.ls.*id2489385430.pu.Omegab.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,0,0;
id2489385430.pu.Omegab.*(-id2489385430.ss.Omegag0+id2489385430.ss.Omegavsm0),-1./id2489385430.parameters.lf.*(id2489385430.parameters.kpc+id2489385430.parameters.rf).*id2489385430.pu.Omegab,0,(-1-id2489385430.parameters.kad+id2489385430.parameters.kffv)./id2489385430.parameters.lf.*id2489385430.pu.Omegab,0,0,0,id2489385430.parameters.kic./id2489385430.parameters.lf.*id2489385430.pu.Omegab,0,id2489385430.parameters.kad./id2489385430.parameters.lf.*id2489385430.pu.Omegab,id2489385430.parameters.kdrpq.*id2489385430.parameters.kpc.*id2489385430.parameters.ls.*id2489385430.pu.Omegab.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),1./id2489385430.parameters.lf.*id2489385430.pu.Omegab.*(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2).^(-2).*(id2489385430.parameters.kdrpq.*id2489385430.parameters.kpc.*id2489385430.parameters.ls.*(id2489385430.ss.qm0-id2489385430.ss.qref0).*(id2489385430.parameters.rs.^2-id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2)+id2489385430.ss.ild0.*id2489385430.parameters.lf.*(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2).^2+id2489385430.parameters.kpc.*id2489385430.parameters.ls.*(id2489385430.parameters.rs.^2.*(id2489385430.ss.vomd0-id2489385430.ss.vref0)+2.*id2489385430.parameters.ls.*id2489385430.parameters.rs.*id2489385430.ss.vomq0.*id2489385430.ss.Omegavsm0+id2489385430.parameters.ls.^2.*(-id2489385430.ss.vomd0+id2489385430.ss.vref0).*id2489385430.ss.Omegavsm0.^2)),0,id2489385430.parameters.kpc.*id2489385430.parameters.ls.*id2489385430.pu.Omegab.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),-id2489385430.parameters.kpc.*id2489385430.parameters.rs.*id2489385430.pu.Omegab./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,0,0;
1./id2489385430.parameters.cf.*id2489385430.pu.Omegab,0,0,id2489385430.pu.Omegab.*id2489385430.ss.Omegag0,-1./id2489385430.parameters.cf.*id2489385430.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,1./id2489385430.parameters.cf.*id2489385430.pu.Omegab,-id2489385430.pu.Omegab.*id2489385430.ss.Omegag0,0,0,-1./id2489385430.parameters.cf.*id2489385430.pu.Omegab,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,1./id2489385430.parameters.lg.*id2489385430.pu.Omegab,0,-1./id2489385430.parameters.lg.*id2489385430.parameters.rg.*id2489385430.pu.Omegab,id2489385430.pu.Omegab.*id2489385430.ss.Omegag0,0,0,0,0,0,0,1./id2489385430.parameters.lg.*(-id2489385430.ss.vnq0.*id2489385430.pu.Omegab.*cos(id2489385430.ss.DeltaThetavsm0)+id2489385430.ss.vnd0.*id2489385430.pu.Omegab.*sin(id2489385430.ss.DeltaThetavsm0)),0,0,0,0,0,0;
0,0,0,1./id2489385430.parameters.lg.*id2489385430.pu.Omegab,-id2489385430.pu.Omegab.*id2489385430.ss.Omegag0,-1./id2489385430.parameters.lg.*id2489385430.parameters.rg.*id2489385430.pu.Omegab,0,0,0,0,0,0,1./id2489385430.parameters.lg.*id2489385430.pu.Omegab.*(id2489385430.ss.vnd0.*cos(id2489385430.ss.DeltaThetavsm0)+id2489385430.ss.vnq0.*sin(id2489385430.ss.DeltaThetavsm0)),0,0,0,0,0,0;
-1,0,0,0,0,0,0,0,0,0,-id2489385430.parameters.kdrpq.*id2489385430.parameters.rs./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),id2489385430.parameters.ls.*(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2).^(-2).*(-id2489385430.parameters.rs.^2.*id2489385430.ss.vomq0+2.*id2489385430.parameters.ls.*id2489385430.parameters.rs.*(id2489385430.parameters.kdrpq.*(id2489385430.ss.qm0-id2489385430.ss.qref0)+id2489385430.ss.vomd0-id2489385430.ss.vref0).*id2489385430.ss.Omegavsm0+id2489385430.parameters.ls.^2.*id2489385430.ss.vomq0.*id2489385430.ss.Omegavsm0.^2),0,-id2489385430.parameters.rs./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),-id2489385430.parameters.ls.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,0,0;
0,-1,0,0,0,0,0,0,0,0,id2489385430.parameters.kdrpq.*id2489385430.parameters.ls.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),id2489385430.parameters.ls.*(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2).^(-2).*(id2489385430.parameters.rs.^2.*(id2489385430.ss.vomd0-id2489385430.ss.vref0)+2.*id2489385430.parameters.ls.*id2489385430.parameters.rs.*id2489385430.ss.vomq0.*id2489385430.ss.Omegavsm0+id2489385430.parameters.ls.^2.*(-id2489385430.ss.vomd0+id2489385430.ss.vref0).*id2489385430.ss.Omegavsm0.^2+id2489385430.parameters.kdrpq.*(id2489385430.ss.qm0-id2489385430.ss.qref0).*(id2489385430.parameters.rs.^2-id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2)),0,id2489385430.parameters.ls.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),-id2489385430.parameters.rs./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,0,0;
0,0,id2489385430.parameters.Omegaad,0,0,0,0,0,-id2489385430.parameters.Omegaad,0,0,0,0,0,0,0,0,0,0;
0,0,0,id2489385430.parameters.Omegaad,0,0,0,0,0,-id2489385430.parameters.Omegaad,0,0,0,0,0,0,0,0,0;
0,0,-id2489385430.ss.ioq0.*id2489385430.parameters.Omegaqf,id2489385430.ss.iod0.*id2489385430.parameters.Omegaqf,id2489385430.ss.voq0.*id2489385430.parameters.Omegaqf,-id2489385430.ss.vod0.*id2489385430.parameters.Omegaqf,0,0,0,0,-id2489385430.parameters.Omegaqf,0,0,0,0,0,0,0,0;
0,0,-id2489385430.ss.iod0./id2489385430.parameters.ta,-id2489385430.ss.ioq0./id2489385430.parameters.ta,-1./id2489385430.parameters.ta.*id2489385430.ss.vod0,-1./id2489385430.parameters.ta.*id2489385430.ss.voq0,0,0,0,0,0,-(id2489385430.parameters.kd+id2489385430.parameters.kdrpOmega)./id2489385430.parameters.ta,0,0,0,-id2489385430.parameters.kd.*id2489385430.parameters.kppll.*id2489385430.ss.vpllq0./(id2489385430.parameters.ta.*id2489385430.ss.vplld0.^2+id2489385430.parameters.ta.*id2489385430.ss.vpllq0.^2),id2489385430.parameters.kd.*id2489385430.parameters.kppll.*id2489385430.ss.vplld0./(id2489385430.parameters.ta.*id2489385430.ss.vplld0.^2+id2489385430.parameters.ta.*id2489385430.ss.vpllq0.^2),id2489385430.parameters.kd.*id2489385430.parameters.kipll./id2489385430.parameters.ta,0;
0,0,0,0,0,0,0,0,0,0,0,id2489385430.pu.Omegab,0,0,0,0,0,0,0;
0,0,id2489385430.parameters.Omegavo,0,0,0,0,0,0,0,0,0,0,-id2489385430.parameters.Omegavo,0,0,0,0,0;
0,0,0,id2489385430.parameters.Omegavo,0,0,0,0,0,0,0,0,0,0,-id2489385430.parameters.Omegavo,0,0,0,0;
0,0,id2489385430.parameters.Omegalppll.*cos(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0),id2489385430.parameters.Omegalppll.*sin(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0),0,0,0,0,0,0,0,0,-id2489385430.ss.voq0.*id2489385430.parameters.Omegalppll.*cos(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0)+id2489385430.ss.vod0.*id2489385430.parameters.Omegalppll.*sin(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0),0,0,-id2489385430.parameters.Omegalppll,0,0,id2489385430.parameters.Omegalppll.*(id2489385430.ss.voq0.*cos(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0)-id2489385430.ss.vod0.*sin(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0));
0,0,-id2489385430.parameters.Omegalppll.*sin(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0),id2489385430.parameters.Omegalppll.*cos(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0),0,0,0,0,0,0,0,0,id2489385430.parameters.Omegalppll.*(id2489385430.ss.vod0.*cos(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0)+id2489385430.ss.voq0.*sin(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0)),0,0,0,-id2489385430.parameters.Omegalppll,0,-id2489385430.parameters.Omegalppll.*(id2489385430.ss.vod0.*cos(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0)+id2489385430.ss.voq0.*sin(id2489385430.ss.DeltaThetapll0-id2489385430.ss.DeltaThetavsm0));
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-id2489385430.ss.vpllq0./(id2489385430.ss.vplld0.^2+id2489385430.ss.vpllq0.^2),id2489385430.ss.vplld0./(id2489385430.ss.vplld0.^2+id2489385430.ss.vpllq0.^2),0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-id2489385430.parameters.kppll.*id2489385430.ss.vpllq0./(id2489385430.ss.vplld0.^2+id2489385430.ss.vpllq0.^2).*id2489385430.pu.Omegab,id2489385430.parameters.kppll.*id2489385430.ss.vplld0./(id2489385430.ss.vplld0.^2+id2489385430.ss.vpllq0.^2).*id2489385430.pu.Omegab,id2489385430.parameters.kipll.*id2489385430.pu.Omegab,0];

id2489385430.B=[0,id2489385430.parameters.kdrpq.*id2489385430.parameters.kpc.*id2489385430.parameters.rs.*id2489385430.pu.Omegab./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,id2489385430.parameters.kpc.*id2489385430.parameters.rs.*id2489385430.pu.Omegab./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,id2489385430.ss.ilq0.*id2489385430.pu.Omegab;
0,-id2489385430.parameters.kdrpq.*id2489385430.parameters.kpc.*id2489385430.parameters.ls.*id2489385430.pu.Omegab.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,-id2489385430.parameters.kpc.*id2489385430.parameters.ls.*id2489385430.pu.Omegab.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.lf.*id2489385430.parameters.rs.^2+id2489385430.parameters.lf.*id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,-id2489385430.ss.ild0.*id2489385430.pu.Omegab;
0,0,0,0,0,0,id2489385430.ss.voq0.*id2489385430.pu.Omegab;
0,0,0,0,0,0,-id2489385430.ss.vod0.*id2489385430.pu.Omegab;
0,0,-1./id2489385430.parameters.lg.*id2489385430.pu.Omegab.*cos(id2489385430.ss.DeltaThetavsm0),-1./id2489385430.parameters.lg.*id2489385430.pu.Omegab.*sin(id2489385430.ss.DeltaThetavsm0),0,0,id2489385430.ss.ioq0.*id2489385430.pu.Omegab;
0,0,1./id2489385430.parameters.lg.*id2489385430.pu.Omegab.*sin(id2489385430.ss.DeltaThetavsm0),-1./id2489385430.parameters.lg.*id2489385430.pu.Omegab.*cos(id2489385430.ss.DeltaThetavsm0),0,0,-id2489385430.ss.iod0.*id2489385430.pu.Omegab;
0,id2489385430.parameters.kdrpq.*id2489385430.parameters.rs./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,id2489385430.parameters.rs./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0;
0,-id2489385430.parameters.kdrpq.*id2489385430.parameters.ls.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0,-id2489385430.parameters.ls.*id2489385430.ss.Omegavsm0./(id2489385430.parameters.rs.^2+id2489385430.parameters.ls.^2.*id2489385430.ss.Omegavsm0.^2),0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
1./id2489385430.parameters.ta,0,0,0,0,id2489385430.parameters.kdrpOmega./id2489385430.parameters.ta,0;
0,0,0,0,0,0,-id2489385430.pu.Omegab;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,-id2489385430.pu.Omegab];

id2489385430.C = eye(size(id2489385430.A,1));
id2489385430.sys = ss(id2489385430.A, id2489385430.B, id2489385430.C, 0, 'statename', {strcat(id1947728959,'_ild') strcat(id1947728959,'_ilq') strcat(id1947728959,'_vod') strcat(id1947728959,'_voq') strcat(id1947728959,'_iod') strcat(id1947728959,'_ioq') strcat(id1947728959,'_Gammad') strcat(id1947728959,'_Gammaq') strcat(id1947728959,'_Phid') strcat(id1947728959,'_Phiq') strcat(id1947728959,'_qm') strcat(id1947728959,'_Omegavsm') strcat(id1947728959,'_DeltaThetavsm') strcat(id1947728959,'_vomd') strcat(id1947728959,'_vomq') strcat(id1947728959,'_vplld') strcat(id1947728959,'_vpllq') strcat(id1947728959,'_Epsilonpll') strcat(id1947728959,'_DeltaThetapll') },'inputname', {strcat(id1947728959,'_pref') strcat(id1947728959,'_qref') strcat(id1947728959,'_vnd') strcat(id1947728959,'_vnq') strcat(id1947728959,'_vref') strcat(id1947728959,'_Omegaref') strcat(id1947728959,'_Omegag') },'outputname', {strcat(id1947728959,'_ild') strcat(id1947728959,'_ilq') strcat(id1947728959,'_vod') strcat(id1947728959,'_voq') strcat(id1947728959,'_iod') strcat(id1947728959,'_ioq') strcat(id1947728959,'_Gammad') strcat(id1947728959,'_Gammaq') strcat(id1947728959,'_Phid') strcat(id1947728959,'_Phiq') strcat(id1947728959,'_qm') strcat(id1947728959,'_Omegavsm') strcat(id1947728959,'_DeltaThetavsm') strcat(id1947728959,'_vomd') strcat(id1947728959,'_vomq') strcat(id1947728959,'_vplld') strcat(id1947728959,'_vpllq') strcat(id1947728959,'_Epsilonpll') strcat(id1947728959,'_DeltaThetapll') });
end

function id3368576760=id3325638676(id3627625384,id2436168768)
id3368576760= [1./id2436168768.parameters.ls.*id2436168768.pu.Omegab.*(-id3627625384(2).*id2436168768.parameters.rs+id2436168768.ss.vref0-id2436168768.parameters.cf.*id2436168768.parameters.rs.*id2436168768.ss.Omegag0./(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).*(id3627625384(2).*id2436168768.ss.Omegag0.*(-id2436168768.parameters.lg.*id2436168768.parameters.rs+id2436168768.parameters.rg.*(id2436168768.parameters.ls+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.parameters.rs)+id2436168768.parameters.cf.*id2436168768.parameters.lg.^2.*id2436168768.parameters.rs.*id2436168768.ss.Omegag0.^2)+id2436168768.parameters.rs.*(id2436168768.parameters.cf.*id2436168768.parameters.rg.* ...
  id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.ss.vnq0.*(-1+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.Omegag0.^2)).*cos(id3627625384(1))+id2436168768.parameters.rs.*(id2436168768.ss.vnd0+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2).*sin(id3627625384(1)))-id2436168768.parameters.ls.*id2436168768.ss.Omegag0./(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+ ...
  id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).*(id3627625384(2).*id2436168768.ss.Omegag0.*(id2436168768.parameters.lg+id2436168768.parameters.ls+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.parameters.rs-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^2)+(id2436168768.ss.vnq0+id2436168768.parameters.cf.*id2436168768.parameters.rs.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0-id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2).*cos(id3627625384(1))+(-id2436168768.ss.vnd0+id2436168768.parameters.cf.*id2436168768.parameters.rs.*id2436168768.ss.vnq0.* ...
  id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2).*sin(id3627625384(1)))+1./(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).*(id3627625384(2).*(-id2436168768.parameters.rg.*(id2436168768.parameters.rg+id2436168768.parameters.rs)-id2436168768.parameters.lg.*(id2436168768.parameters.lg+id2436168768.parameters.ls).*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.^2.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.* ...
  id2436168768.parameters.lg.^2.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^4)+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).*id2436168768.ss.vnd0-id2436168768.parameters.lg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^3).*cos(id3627625384(1))+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).*id2436168768.ss.vnq0+id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.* ...
  id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^3).*sin(id3627625384(1)))-id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^2./(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).*(id3627625384(2).*(-id2436168768.parameters.rg.*(id2436168768.parameters.rg+id2436168768.parameters.rs)-id2436168768.parameters.lg.*(id2436168768.parameters.lg+id2436168768.parameters.ls).*id2436168768.ss.Omegag0.^2+ ...
  id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.^2.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.lg.^2.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^4)+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).*id2436168768.ss.vnd0-id2436168768.parameters.lg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^3).*cos(id3627625384(1))+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).* ...
  id2436168768.ss.vnq0+id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^3).*sin(id3627625384(1)))+id2436168768.parameters.kdrpq.*(id2436168768.ss.qref0+id3627625384(2)./(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).*(id3627625384(2).* ...
  id2436168768.ss.Omegag0.*(-id2436168768.parameters.lg.*id2436168768.parameters.rs+id2436168768.parameters.rg.*(id2436168768.parameters.ls+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.parameters.rs)+id2436168768.parameters.cf.*id2436168768.parameters.lg.^2.*id2436168768.parameters.rs.*id2436168768.ss.Omegag0.^2)+id2436168768.parameters.rs.*(id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.ss.vnq0.*(-1+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.Omegag0.^2)).*cos(id3627625384(1))+id2436168768.parameters.rs.*(id2436168768.ss.vnd0+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0-id2436168768.parameters.cf.* ...
  id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2).*sin(id3627625384(1)))+(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).^(-2).*(id3627625384(2).*id2436168768.ss.Omegag0.*(id2436168768.parameters.lg+id2436168768.parameters.ls+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.parameters.rs-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^2)+(id2436168768.ss.vnq0+id2436168768.parameters.cf.*id2436168768.parameters.rs.*id2436168768.ss.vnd0.* ...
  id2436168768.ss.Omegag0-id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2).*cos(id3627625384(1))+(-id2436168768.ss.vnd0+id2436168768.parameters.cf.*id2436168768.parameters.rs.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2).*sin(id3627625384(1))).*(id3627625384(2).*(-id2436168768.parameters.rg.*(id2436168768.parameters.rg+id2436168768.parameters.rs)-id2436168768.parameters.lg.*(id2436168768.parameters.lg+id2436168768.parameters.ls).*id2436168768.ss.Omegag0.^2+ ...
  id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.^2.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.lg.^2.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^4)+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).*id2436168768.ss.vnd0-id2436168768.parameters.lg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^3).*cos(id3627625384(1))+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).* ...
  id2436168768.ss.vnq0+id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^3).*sin(id3627625384(1))))),1./id2436168768.parameters.ta.*(id2436168768.ss.pref0-id2436168768.parameters.kdrpOmega.*id2436168768.ss.Omegag0+id2436168768.parameters.kdrpOmega.*id2436168768.ss.Omegaref0-(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+ ...
  id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).^(-2).*(id3627625384(2).*id2436168768.ss.Omegag0.*(-id2436168768.parameters.lg.*id2436168768.parameters.rs+id2436168768.parameters.rg.*(id2436168768.parameters.ls+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.parameters.rs)+id2436168768.parameters.cf.*id2436168768.parameters.lg.^2.*id2436168768.parameters.rs.*id2436168768.ss.Omegag0.^2)+id2436168768.parameters.rs.*(id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.ss.vnq0.*(-1+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.Omegag0.^2)).*cos(id3627625384(1))+ ...
  id2436168768.parameters.rs.*(id2436168768.ss.vnd0+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2).*sin(id3627625384(1))).*(id3627625384(2).*id2436168768.ss.Omegag0.*(id2436168768.parameters.lg+id2436168768.parameters.ls+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.parameters.rs-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^2)+(id2436168768.ss.vnq0+id2436168768.parameters.cf.*id2436168768.parameters.rs.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0-id2436168768.parameters.cf.*id2436168768.parameters.ls.* ...
  id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2).*cos(id3627625384(1))+(-id2436168768.ss.vnd0+id2436168768.parameters.cf.*id2436168768.parameters.rs.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2).*sin(id3627625384(1)))+id3627625384(2)./(id2436168768.parameters.rg+id2436168768.parameters.rs-id2436168768.parameters.cf.*(id2436168768.parameters.ls.*id2436168768.parameters.rg+id2436168768.parameters.lg.*id2436168768.parameters.rs).*id2436168768.ss.Omegag0.^2).*(id3627625384(2).*(-id2436168768.parameters.rg.*(id2436168768.parameters.rg+ ...
  id2436168768.parameters.rs)-id2436168768.parameters.lg.*(id2436168768.parameters.lg+id2436168768.parameters.ls).*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.^2.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.lg.^2.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^4)+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).*id2436168768.ss.vnd0-id2436168768.parameters.lg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^3).*cos( ...
  id3627625384(1))+(-(id2436168768.parameters.rg+id2436168768.parameters.rs).*id2436168768.ss.vnq0+id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^3).*sin(id3627625384(1)))),1./(id2436168768.parameters.rs.*(-1+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.Omegag0.^2)+id2436168768.parameters.rg.*(-1+id2436168768.parameters.cf.*id2436168768.parameters.ls.* ...
  id2436168768.ss.Omegag0.^2)).*(cos(id3627625384(3)-id3627625384(1)).*(id3627625384(2).*id2436168768.ss.Omegag0.*(id2436168768.parameters.ls.*id2436168768.parameters.rg-id2436168768.parameters.lg.*id2436168768.parameters.rs+id2436168768.parameters.cf.*id2436168768.parameters.rg.^2.*id2436168768.parameters.rs+id2436168768.parameters.cf.*id2436168768.parameters.lg.^2.*id2436168768.parameters.rs.*id2436168768.ss.Omegag0.^2)+id2436168768.parameters.rs.*(id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0+id2436168768.ss.vnq0.*(-1+id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.Omegag0.^2)) ...
  .*cos(id3627625384(1))+id2436168768.parameters.rs.*(id2436168768.ss.vnd0+id2436168768.parameters.cf.*id2436168768.parameters.rg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2).*sin(id3627625384(1)))+sin(id3627625384(3)-id3627625384(1)).*(id3627625384(2).*(id2436168768.parameters.rg.*id2436168768.parameters.rs+id2436168768.parameters.rg.^2.*(1-id2436168768.parameters.cf.*id2436168768.parameters.ls.* ...
  id2436168768.ss.Omegag0.^2)+id2436168768.parameters.lg.*id2436168768.ss.Omegag0.^2.*(id2436168768.parameters.lg+id2436168768.parameters.ls-id2436168768.parameters.cf.*id2436168768.parameters.lg.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^2))+(id2436168768.parameters.rs.*id2436168768.ss.vnd0+id2436168768.parameters.lg.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.*(1-id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^2)+id2436168768.parameters.rg.*(id2436168768.ss.vnd0-id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.^2)).*cos(id3627625384(1))+(id2436168768.parameters.rs.*id2436168768.ss.vnq0+ ...
  id2436168768.parameters.lg.*id2436168768.ss.vnd0.*id2436168768.ss.Omegag0.*(-1+id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.Omegag0.^2)+id2436168768.parameters.rg.*(id2436168768.ss.vnq0-id2436168768.parameters.cf.*id2436168768.parameters.ls.*id2436168768.ss.vnq0.*id2436168768.ss.Omegag0.^2)).*sin(id3627625384(1))))];
end