
clc; clear all; close all

h_cut = 1.055e-34;
eV = 1.602e-19;
Armst=1e-10;


t_eV=2.7;
a0_A0=1.42;

t=t_eV*eV;
a0=a0_A0*Armst;
a=3*a0/2;
d12=2*a;

na=30;
n=na-2;
d=6+(n-1)*2;
sn=[1:6;[2:6,1]];
tn=[5;2];
for m=2:n
    nm=(7+(m-2)*2):(8+(m-2)*2);
    if m==2
        c1=4;
    else
        c1=nm(1)-1;
    end
    nc=[nm(1:end-1),c1;nm(2:end),nm(1)];
    sn=[sn,nc];
    if m==2
        tn=[tn,[nm(end);3]];
    else
        tn=[tn,[nm(end);nm(1)-2]];
    end
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
xlabel('k_x \rightarrow');ylabel('E(k_x) in eV \rightarrow')
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi/d','-\pi/2d','0','\pi/2d','\pi/d'})
title(['AGNR : Energy versus k_x fro N_a=',num2str(na)])



res=1e3;
kv=linspace(-pi/a,pi/a,res);
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
xlabel('k_x in radian/m \rightarrow');ylabel('E(k_x) in eV \rightarrow')
title(['AGNR : Energy versus k_x fro N_a=',num2str(na)])


