function G=Elastic_GreenTensor_Dhalf_SIP(omega,kp,ks,y,x)
% ������ʽ����
%y is point source
tic
% 2016 12 16 
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

% ���ü�ӱ���

mus=@(t) (ks^2-t.^2).^(1/2);
mup=@(t) (kp^2-t.^2).^(1/2);
delta=@(t) t.^2+mus(t).*mup(t);
forec=@(t) 1i*w2./(2*pi*delta(t));

Ess=@(t) exp(1i*mus(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
Epp=@(t) exp(1i*mup(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
Esp=@(t) exp(1i*(mus(t)*x(2,:)+mup(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));
Eps=@(t) exp(1i*(mup(t)*x(2,:)+mus(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));
B1=@(t) t.*t.*mus(t);
B2=@(t) -t.^3;
B3=@(t) -t.*mus(t).*mup(t);
B4=@(t) t.*t.*mup(t);
A1=@(t) forec(t).*(B1(t).*Ess(t)+B1(t).*Epp(t)-B1(t).*Esp(t)-B1(t).*Eps(t));
A2=@(t) forec(t).*(B2(t).*Ess(t)-B3(t).*Epp(t)-B2(t).*Esp(t)+B3(t).*Eps(t));
A3=@(t) forec(t).*(B3(t).*Ess(t)-B2(t).*Epp(t)-B3(t).*Esp(t)+B2(t).*Eps(t));
A4=@(t) forec(t).*(B4(t).*Ess(t)+B4(t).*Epp(t)-B4(t).*Esp(t)-B4(t).*Eps(t));

Gc1=integral(@(t) A1(t) ,-l,l, 'ArrayValued',true);
Gc2=integral(@(t) A2(t) ,-l,l, 'ArrayValued',true);
Gc3=integral(@(t) A3(t) ,-l,l, 'ArrayValued',true);
Gc4=integral(@(t) A4(t) ,-l,l, 'ArrayValued',true);

G(1,:)=G1(1,:)-G2(1,:)+Gc1;
G(2,:)=G1(2,:)-G2(2,:)+Gc2;
G(3,:)=G1(3,:)-G2(3,:)+Gc3;
G(4,:)=G1(4,:)-G2(4,:)+Gc4;
toc