clc; clear all; close all

h_cut = 1.055e-34;
eV = 1.602e-19;
Armst=1e-10;

ep0_eV=0;
t_eV=2.7;
a0_A0=1.42;

ep0=ep0_eV*eV;
t=t_eV*eV;
a0=a0_A0*Armst;

a=3*a0/2;
b=sqrt(3)*a0/2;
a1=[a;b];
a2=[a;-b];

res=1e3;
rang=2*pi/a;
kx=linspace(-rang,rang,res);
ky=zeros(1,res);
kv=[kx; ky];
Ev=zeros(size(kv));

for i=1:length(kv)
    k=kv(:,i);
    h0=-t*(1+2*exp(1i*k(1)*a)*cos(k(2)*b));
    h0a=abs(h0);
    Ev(:,i)=[-h0a;h0a];
end

figure
plot(kx*a,Ev(1,:)/eV, 'LineWidth', 2)
grid on;hold on
plot(kx*a,Ev(2,:)/eV, 'LineWidth', 2)
xlabel('k_x \rightarrow');ylabel('E(\bf{k}) in eV \rightarrow')
xticks([-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'-2\pi/a','-3\pi/2a','-\pi/a','-\pi/2a','0','\pi/2a','\pi/a','3\pi/2a', '2\pi/a'})
title('Energy versus k_x ; for k_y=0')
legend('Valence band','Conduction band');

rang=2*pi/b;
kx=zeros(1,res);
ky=linspace(-rang,rang,res);
kv=[kx; ky];
Ev=zeros(size(kv));

for i=1:res
    k=kv(:,i);
    h0=-t*(1+2*exp(1i*k(1)*a)*cos(k(2)*b));
    h0a=abs(h0);
    Ev(:,i)=[-h0a;h0a];
end

figure
plot(ky*b,Ev(1,:)/eV, 'LineWidth', 2)
grid on;hold on
plot(ky*b,Ev(2,:)/eV, 'LineWidth', 2)
xlabel('k_y \rightarrow');ylabel('E(\bf{k}) in eV \rightarrow')
xticks([-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'-2\pi/b','-3\pi/2b','-\pi/b','-\pi/2b','0','\pi/2b','\pi/b','3\pi/2b', '2\pi/b'})
title('Energy versus k_y ; for k_x=0')
legend('Valence band','Conduction band');


res=2e2;
kx=linspace(-6/5*pi/a,6/5*pi/a,res);
ky=linspace(-6/5*pi/b,6/5*pi/b,res);
[Kx,Ky]=meshgrid(kx,ky);
kv=[Kx(:).';Ky(:).'];
Ev=zeros(size(kv));

for i=1:length(kv)
    k=kv(:,i);
    h0=-t*(1+2*exp(1i*k(1)*a)*cos(k(2)*b));
    h0a=abs(h0);
    Ev(:,i)=[-h0a;h0a];
end

E1=reshape(Ev(1,:),res,res);
E2=reshape(Ev(2,:),res,res);

figure
h1=surf(Kx*a,Ky*b,E1/eV);
hold on
h2=surf(Kx*a,Ky*b,E2/eV);
shading interp
set(h1,'edgecolor','k');set(h2,'edgecolor','k')
xlabel('k_x \rightarrow');ylabel('k_y \rightarrow');zlabel('E(\bf{k}) in eV \rightarrow')
xticks([-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
yticks([-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'-2\pi/a','-3\pi/2a','-\pi/a','-\pi/2a','0','\pi/2a','\pi/a','3\pi/2a', '2\pi/a'})
yticklabels({'-2\pi/b','-3\pi/2b','-\pi/b','-\pi/2b','0','\pi/2b','\pi/b','3\pi/2b', '2\pi/b'})
title('E-k surface plot')
legend('Valence band','Conduction band');

figure
pcolor(Kx*a,Ky*b,E2/eV);
shading('interp');
hold on
contour(Kx*a,Ky*b,E2/eV, 'LineColor', 'b');
colorbar
xlabel('k_x \rightarrow');ylabel('k_y \rightarrow');zlabel('E(\bf{k}) in eV \rightarrow')
xticks([-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
yticks([-2*pi -3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'-2\pi/a','-3\pi/2a','-\pi/a','-\pi/2a','0','\pi/2a','\pi/a','3\pi/2a', '2\pi/a'})
yticklabels({'-2\pi/b','-3\pi/2b','-\pi/b','-\pi/2b','0','\pi/2b','\pi/b','3\pi/2b', '2\pi/b'})
title('E-k color and contour plot of conduction band')


res=60;
kx=linspace(-pi/a,pi/a,res);
ky=linspace(-pi/b,pi/b,res);
[Kx,Ky]=meshgrid(kx,ky);
kv=[Kx(:).';Ky(:).'];
Ev=zeros(size(kv));

for i=1:length(kv)
    k=kv(:,i);
    h0=-t*(1+2*exp(1i*k(1)*a)*cos(k(2)*b));
    h0a=abs(h0);
    Ev(:,i)=[-h0a;h0a];
end

E1=reshape(Ev(1,:),res,res);
E2=reshape(Ev(2,:),res,res);

figure
h1=surf(Kx,Ky,E1);
hold on
h2=surf(Kx,Ky,E2);
shading interp
set(h1,'edgecolor','k');set(h2,'edgecolor','k')
xlabel('k_x \rightarrow');ylabel('k_y \rightarrow');zlabel('E(\bf{k}) in eV \rightarrow')
title('E-k surface plot for bulk graphene')