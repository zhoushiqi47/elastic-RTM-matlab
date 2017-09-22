n=100;
b=linspace(0,19,100);
t=pi/2-1i*b;
x=20;
ks=2;
kp=1;
phi=pi/4*ones(1,100);
(sin(t+phi).^2)./(sin(t+phi).^2+cos(t+phi).*(ks^2-sin(t+phi).^2).^(1/2)).*cos(t+phi).*exp(1i*20.*COSC(t))