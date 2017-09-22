function G=Elastic_GreenTensor_Thalf_SIP_test(omega,kp,ks,y,x)
% 向量形式积分
%y is point source
% 2016 9 21 
w2 = 1/omega^2 ;
n=size(x);
G=zeros(4,n(2));

l=10*(kp+ks);
yy=y;
yy(2,:)=-y(2,:);
G1=Elastic_GreenTensor_2D(omega,kp,ks,y,x);
G2=Elastic_GreenTensor_2D(omega,kp,ks,yy,x);

kp=kp*(1+0.0000001i);
ks=ks*(1+0.0000001i);

% 设置间接变量

mus=@(t) (ks^2-t.^2).^(1/2);
mup=@(t) (kp^2-t.^2).^(1/2);
beta=@(t) ks^2-2.*t.^2;
delta=@(t) beta(t).^2 + 4*t.^2.*mus(t).*mup(t);
forec=@(t) 1i*w2./(2*pi*delta(t));

Ess=@(t) exp(1i*mus(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
Epp=@(t) exp(1i*mup(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
A=@(t) [t;t];
AA=@(t) A(t)*t;
B=@(t) [t,t];
C=@(t) B(t)*AA(t);

A1=@(t) forec(t).*(2*ks*ks*C(t).*mus(t)*mup(t).*t.*Ess(t)-     ks*ks.*beta(t).*t.*Epp(t));


Gc1=integral(@(t) A1(t) ,-l,l, 'ArrayValued',true);


G=Gc1;

