function [Z_frd]=IM_C(f,f1,C)

omega=f*2*pi;
n=1;
for s=1i*omega
    Zp=1/(C*(s+f1*2*pi*1i));
    Zn=1/(C*(s-f1*2*pi*1i));
    Z(:,:,n)=[Zp 0; 0 Zn];
    n=n+1;
end

%%
Z_frd=frd(Z,f,'FrequencyUnit','Hz');