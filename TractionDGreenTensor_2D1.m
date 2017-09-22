function G=TractionDGreenTensor_2D1(omega,kp,ks,y,x)
%% n=(0,1)
%% operator T to x
%% compute the traction green tensor of elastic wave operator in 2D
%% shear wave velocity
%% x2=0
%% pressure wave Green Tensor


tic

% 向量形式积分
%y is point source
% 2016 12 16 

n=size(x);
G=zeros(4,n(2));

l=10*(kp+ks);
kp=kp*(1+0.0000001i);
ks=ks*(1+0.0000001i);

% 设置间接变量

mus=@(t) (ks^2-t.^2).^(1/2);
mup=@(t) (kp^2-t.^2).^(1/2);
delta=@(t) t.^2+mus(t).*mup(t);
forec=@(t) 1./(2*pi*delta(t));

Ess=@(t) exp(1i*mus(t)*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Epp=@(t) exp(1i*mup(t)*(y(2,:))+1i*t*(x(1,:)-y(1,:)));


B1=@(t) mus(t).*mup(t);
B2=@(t) t.*mus(t);
B3=@(t) t.*mup(t);
B4=@(t) t.*t;



A1=@(t) forec(t).*(B1(t).*Ess(t)+B4(t).*Epp(t));
A2=@(t) forec(t).*(B2(t).*Ess(t)-B2(t).*Epp(t));
A3=@(t) forec(t).*(B3(t).*Ess(t)-B3(t).*Epp(t));
A4=@(t) forec(t).*(B4(t).*Ess(t)+B1(t).*Epp(t));

Gc1=integral(@(t) A1(t) ,-l,l, 'ArrayValued',true);
Gc2=integral(@(t) A2(t) ,-l,l, 'ArrayValued',true);
Gc3=integral(@(t) A3(t) ,-l,l, 'ArrayValued',true);
Gc4=integral(@(t) A4(t) ,-l,l, 'ArrayValued',true);

G(1,:)=Gc1;
G(2,:)=Gc2;
G(3,:)=Gc3;
G(4,:)=Gc4;
toc