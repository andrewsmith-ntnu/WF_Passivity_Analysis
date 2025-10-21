function [Z_frd]=IM_RL(f,f1,R,L)
%Function to calculate the sequence-domain impedance of an RL element

%   f = array of frequencies, in Hz
%   f1 = fundamental frequency, in Hz
%   R = resistance
%   L = inductance

omega=f*2*pi;
n=1;
for s=1i*omega
    Zp=R+L*(s+f1*2*pi*1i);
    Zn=R+L*(s-f1*2*pi*1i);
    Z(:,:,n)=[Zp 0; 0 Zn];
    n=n+1;
end

%%
Z_frd=frd(Z,f,'FrequencyUnit','Hz');