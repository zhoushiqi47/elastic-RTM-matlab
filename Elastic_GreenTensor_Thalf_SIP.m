function G=Elastic_GreenTensor_Thalf_SIP(omega,kp,ks,y,x)
% 向量形式积分
%y is point source
% 2016 9 21 
tic
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
Esp=@(t) exp(1i*(mus(t)*x(2,:)+mup(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));
Eps=@(t) exp(1i*(mup(t)*x(2,:)+mus(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));

A1=@(t) forec(t).*((mus(t).*beta(t).^2).*Ess(t)+     (4*mus(t).*t.^4).*Epp(t)+        (2*mus(t).*beta(t).*t.^2).*Esp(t)+      (2*mus(t).*beta(t).*t.^2).*Eps(t));
A2=@(t) forec(t).*((-beta(t).^2.*t).*Ess(t)+         (4*mus(t).*mup(t).*t.^3).*Epp(t)+(-2*beta(t).*t.^3).*Esp(t)+             (2*mus(t).*mup(t).*beta(t).*t).*Eps(t));
A3=@(t) forec(t).*((-4*mus(t).*mup(t).*t.^3).*Ess(t)+(beta(t).^2.*t).*Epp(t)+         (-2*mus(t).*mup(t).*beta(t).*t).*Esp(t)+(2*beta(t).*t.^3).*Eps(t));
A4=@(t) forec(t).*((4*mup(t).*t.^4).*Ess(t)+         (mup(t).*beta(t).^2).*Epp(t)+    (2*mup(t).*beta(t).*t.^2).*Esp(t)+      (2*mup(t).*beta(t).*t.^2).*Eps(t));

Gc1=integral(@(t) A1(t) ,-l,l, 'ArrayValued',true);
Gc2=integral(@(t) A2(t) ,-l,l, 'ArrayValued',true);
Gc3=integral(@(t) A3(t) ,-l,l, 'ArrayValued',true);
Gc4=integral(@(t) A4(t) ,-l,l, 'ArrayValued',true);

G(1,:)=G1(1,:)-G2(1,:)+Gc1;
G(2,:)=G1(2,:)-G2(2,:)+Gc2;
G(3,:)=G1(3,:)-G2(3,:)+Gc3;
G(4,:)=G1(4,:)-G2(4,:)+Gc4;
toc

    

