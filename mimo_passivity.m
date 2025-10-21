function pass = mimo_passivity(Z)

Z_dd=Z.ResponseData(1,1,:);
Z_qd=Z.ResponseData(2,1,:);
Z_dq=Z.ResponseData(1,2,:);
Z_qq=Z.ResponseData(2,2,:);

a=Z_dd+conj(Z_dd);
b=Z_qq+conj(Z_qq);
c=Z_dq+conj(Z_qd);

PD=Z;
PD.ResponseData=[a conj(c); c b];

for fr=1:length(PD.Frequency)
    PD_f=PD.ResponseData(:,:,fr);
    sigma=eig(PD_f);
    eig_f(:,fr)=sigma;
end

pass=min(eig_f,[],1);