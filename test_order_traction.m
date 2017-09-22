n=100;
b=10;
x=linspace(200,1000,n);
ks=0.51;
kp=1;
phi=pi/4;
%phi=linspace(pi/8,pi/4,n);
E=@(t) (sin(t+phi).^2)./(sin(t+phi).^2+cos(t+phi).*(ks-sin(t+phi).^2).^(1/2)).*cos(t+phi).*exp(1i*x.*cos(t));
%principle term
y0=(sin(phi).^2)./(sin(phi).^2+cos(phi).*(ks-sin(phi).^2).^(1/2)).*cos(phi).*exp(1i*x).*(-(2*pi*1i)./x).^(1/2);
y=integral(@(t) E(t) ,-pi/4,pi/4, 'ArrayValued',true);

%y1=integral(@(t) E(t) ,-pi/2+b*1i,-pi/2, 'ArrayValued',true);

y2=y-y0;
%y2=y2./x./cos(phi)
z1=log(x);
z2=log(abs(y2));
figure;
plot(z1,z2);
t=(z2(n)-z2(1))/(z1(n)-z1(1))
