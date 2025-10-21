function [cvtr] = NR_steadystate_VSC_GF(cvtr)
%% 'solveNewtonRaphson' function parameters
cnv = sqrt(eps); %Accuracy
Ns = 100;         %Number of iterations
bta = [0.1; 0.2; 0.5; ones(Ns-3,1)]; %First steps
litmax = 20;

cvtr.ss.iod0 = cvtr.ss.pref0;
cvtr.ss.DeltaThetapll0 = -0.01;
cvtr = initialize_NR_VSC_GF(cvtr);

[xsol,dxsol] = solveNewtonRaphson(@(x) function_NR_VSC_GF(x,cvtr),0,cvtr.x0.',cnv,Ns,bta,litmax);
       
cvtr.ss.iod0 = xsol(5);
cvtr.ss.DeltaThetapll0 = xsol(16);
cvtr = initialize_NR_VSC_GF(cvtr);





% Setting three phase initial conditions
cvtr.ss.io0a=(cos(cvtr.ss.DeltaThetapll0)*real(cvtr.ss.io0)-sin(cvtr.ss.DeltaThetapll0)*imag(cvtr.ss.io0))*cvtr.pu.Ib;
cvtr.ss.io0b=(cos(cvtr.ss.DeltaThetapll0-2/3*pi)*real(cvtr.ss.io0)-sin(cvtr.ss.DeltaThetapll0-2/3*pi)*imag(cvtr.ss.io0))*cvtr.pu.Ib;
cvtr.ss.io0c=(cos(cvtr.ss.DeltaThetapll0+2/3*pi)*real(cvtr.ss.io0)-sin(cvtr.ss.DeltaThetapll0+2/3*pi)*imag(cvtr.ss.io0))*cvtr.pu.Ib;

cvtr.ss.il0a=(cos(cvtr.ss.DeltaThetapll0)*real(cvtr.ss.il0)-sin(cvtr.ss.DeltaThetapll0)*imag(cvtr.ss.il0))*cvtr.pu.Ib;
cvtr.ss.il0b=(cos(cvtr.ss.DeltaThetapll0-2/3*pi)*real(cvtr.ss.il0)-sin(cvtr.ss.DeltaThetapll0-2/3*pi)*imag(cvtr.ss.il0))*cvtr.pu.Ib;
cvtr.ss.il0c=(cos(cvtr.ss.DeltaThetapll0+2/3*pi)*real(cvtr.ss.il0)-sin(cvtr.ss.DeltaThetapll0+2/3*pi)*imag(cvtr.ss.il0))*cvtr.pu.Ib;

cvtr.ss.vo0a=(cos(cvtr.ss.DeltaThetapll0)*real(cvtr.ss.vo0)-sin(cvtr.ss.DeltaThetapll0)*imag(cvtr.ss.vo0))*cvtr.pu.Vb;
cvtr.ss.vo0b=(cos(cvtr.ss.DeltaThetapll0-2/3*pi)*real(cvtr.ss.vo0)-sin(cvtr.ss.DeltaThetapll0-2/3*pi)*imag(cvtr.ss.vo0))*cvtr.pu.Vb;
cvtr.ss.vo0c=(cos(cvtr.ss.DeltaThetapll0+2/3*pi)*real(cvtr.ss.vo0)-sin(cvtr.ss.DeltaThetapll0+2/3*pi)*imag(cvtr.ss.vo0))*cvtr.pu.Vb;
