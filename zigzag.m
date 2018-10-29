clc; clear all; close all

h_cut = 1.055e-34;
eV = 1.602e-19;
Armst=1e-10;

t_eV=2.7;
a0_A0=1.42;

t=t_eV*eV;
a0=a0_A0*Armst;
b=sqrt(3)*a0/2;
d12=2*b;

nz=30;
n=nz-1;
d=4+(n-1)*2;
sn=[1:3;2:4];
tn=[1,4;2,3];
for m=2:n
    nm=(5+(m-2)*2):(6+(m-2)*2);
    c1=nm(1)-1;
    nc=[nm(1:end-1),c1;nm(2:end),nm(1)];
    sn=[sn,nc];
    tn=[tn,[nm(end);nm(1)]];
end


res=1e3;
kv=linspace(-pi/d12,pi/d12,res);
Ev=zeros(d,length(kv));

for m=1:length(kv)
    k=kv(m);
    i=[sn(1,:),sn(2,:),tn(1,:),tn(2,:)];
    j=[sn(2,:),sn(1,:),tn(2,:),tn(1,:)];
    v=[repmat(-t, 1,2*size(sn,2)),repmat(-t*exp(1i*k*d12), 1,size(tn,2)),repmat(-t*exp(-1i*k*d12),1, size(tn,2))];
    
    hk=sparse(i,j,v,d,d);
    egv=eig(full(hk));
    Ev(:,m)=egv;
end

figure
plot(kv*d12,Ev/eV, 'b', 'LineWidth',1.5)
xlabel('k_y \rightarrow');ylabel('E(k_y) in eV \rightarrow')
xticks([-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'-2\pi/d','-3\pi/2d','-\pi/d','-\pi/2d','0','\pi/2d','\pi/d','3\pi/2d', '2\pi/d'})
title(['ZGNR : Energy versus k_y fro N_z=',num2str(nz)])



res=1e3;
kv=linspace(-pi/b,pi/b,res);
Ev=zeros(d,length(kv));

for m=1:length(kv)
    k=kv(m);
    i=[sn(1,:),sn(2,:),tn(1,:),tn(2,:)];
    j=[sn(2,:),sn(1,:),tn(2,:),tn(1,:)];
    v=[repmat(-t, 1,2*size(sn,2)),repmat(-t*exp(1i*k*d12), 1,size(tn,2)),repmat(-t*exp(-1i*k*d12),1, size(tn,2))];
    
    hk=sparse(i,j,v,d,d);
    egv=eig(full(hk));
    Ev(:,m)=egv;
end

figure
plot(kv,Ev/eV, 'b', 'LineWidth',1.5)
xlabel('k_y in radian/m \rightarrow');ylabel('E(k_y) in eV \rightarrow')
title(['ZGNR : Energy versus k_y fro N_z=',num2str(nz)])

