kp=pi;
ks=2*pi;
omega=pi;
y=[10;0];
x=[0;10];
yy=[0;-100];
G=Elastic_GreenTensor_Thalf_SIP(omega,kp,ks,y,x)
G1=Elastic_GreenTensor_Thalf_SIP3(omega,kp,ks,y,x)