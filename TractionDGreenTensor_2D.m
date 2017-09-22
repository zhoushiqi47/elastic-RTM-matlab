function G=TractionDGreenTensor_2D(omega,kp,ks,y,x)
%% n=(0,1)
%% operator T to x
%% compute the traction green tensor of elastic wave operator in 2D
%% shear wave velocity
%% pressure wave Green Tensor
tic
w2 = 1/omega^2 ;
mu = omega.^2./(ks.^2);
% 向量形式积分
%y is point source
% 2016 12 16 
n=size(x);
G=zeros(4,n(2));

l=10*(kp+ks);
yy=y;
yy(2,:)=-y(2,:);
G1=TractionGreenTensor_2D(omega,kp,ks,y,x);
G2=TractionGreenTensor_2D(omega,kp,ks,yy,x);

kp=kp*(1+0.0000001i);
ks=ks*(1+0.0000001i);

% 设置间接变量

mus=@(t) (ks^2-t.^2).^(1/2);
mup=@(t) (kp^2-t.^2).^(1/2);
delta=@(t) t.^2+mus(t).*mup(t);
beta=@(t) ks^2-2.*t.^2;
forec=@(t) mu*w2./(2*pi*delta(t));

Ess=@(t) exp(1i*mus(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
Epp=@(t) exp(1i*mup(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
Esp=@(t) exp(1i*(mus(t)*x(2,:)+mup(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));
Eps=@(t) exp(1i*(mup(t)*x(2,:)+mus(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));

B1=@(t) -t.*t.*beta(t);
B2=@(t) 2*t.^3.*mus(t);
B3=@(t) t.*mup(t).*beta(t);
B4=@(t) -2*t.^2.*mus(t).*mup(t);

C2=@(t) -t.*mus(t).*beta(t);
C3=@(t) -2*t.^3.*mup(t);


A1=@(t) forec(t).*(B1(t).*Ess(t)+B4(t).*Epp(t)-B1(t).*Esp(t)-B4(t).*Eps(t));
A2=@(t) forec(t).*(B2(t).*Ess(t)+C2(t).*Epp(t)-B2(t).*Esp(t)-C2(t).*Eps(t));
A3=@(t) forec(t).*(B3(t).*Ess(t)+C3(t).*Epp(t)-B3(t).*Esp(t)-C3(t).*Eps(t));
A4=@(t) forec(t).*(B4(t).*Ess(t)+B1(t).*Epp(t)-B4(t).*Esp(t)-B1(t).*Eps(t));

Gc1=integral(@(t) A1(t) ,-l,l, 'ArrayValued',true);
Gc2=integral(@(t) A2(t) ,-l,l, 'ArrayValued',true);
Gc3=integral(@(t) A3(t) ,-l,l, 'ArrayValued',true);
Gc4=integral(@(t) A4(t) ,-l,l, 'ArrayValued',true);
G(1,:)=G1(1,:)-G2(1,:)+Gc1;
G(2,:)=G1(2,:)-G2(2,:)+Gc2;
G(3,:)=G1(3,:)-G2(3,:)+Gc3;
G(4,:)=G1(4,:)-G2(4,:)+Gc4;
toc