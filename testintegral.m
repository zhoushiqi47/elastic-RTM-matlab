l=1000;
x=zeros(2,l);
y=zeros(2,l);
x(1,:)=20*rand(1,l);
y(1,:)=30*rand(1,l);
x(2,:)=linspace(-1,-10,l);
y(2,:)=linspace(-5,-40,l);
lamda=1/2;
mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

G=Elastic_GreenTensor_Thalf1(omega,kp,ks,x,y);