function G=Elastic_GreenTensor_half2(omega,kp,ks,x,y)
%x is point source%
% 2016 6 22 

w2 = 1/omega^2 ;
k1=kp*(1+0.0000001i);
k2=ks*(1+0.0000001i);


n=size(x);
G=zeros(4,n(2));

l=10*(kp+ks);
xx=x;
xx(2,:)=-x(2,:);
G1=Elastic_GreenTensor_2D(omega,kp,ks,x,y);
G2=Elastic_GreenTensor_2D(omega,kp,ks,xx,y);

mup=@(t) (k1*k1-t.*t).^(1/2);
mus=@(t) (k2*k2-t.*t).^(1/2);
beta=@(t) k2*k2-2*t.*t;
delta=@(t) beta(t).*beta(t)+4*t.*t.*mus(t).*mup(t);
forec=@(t) 1i*w2./(2*pi*delta(t));

ys1=@(t) mus(t).*beta(t).*beta(t).*       exp(1i*(-mus(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
ys2=@(t) t.*beta(t).*beta(t).*            exp(1i*(-mus(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
ys3=@(t) 4*mus(t).*mup(t).*t.^3 .*        exp(1i*(-mus(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
ys4=@(t) 4*mup(t).*t.^4 .*                exp(1i*(-mus(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));

yp1=@(t) 4*mus(t).*t.^4 .*                exp(1i*(-mup(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
yp2=@(t) -4*mus(t).*mup(t).*t.^3 .*       exp(1i*(-mup(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
yp3=@(t) -t.*beta(t).*beta(t) .*          exp(1i*(-mup(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
yp4=@(t) mup(t).*beta(t).*beta(t).*       exp(1i*(-mup(t)*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));

ysp1=@(t) 2*mus(t).*beta(t).*t.^2 .*      exp(1i*(-mup(t)*x(2,:)-mus(t)*y(2,:)+t*(y(1,:)-x(1,:))));
ysp2=@(t) 2*beta(t).*t.^3 .*              exp(1i*(-mup(t)*x(2,:)-mus(t)*y(2,:)+t*(y(1,:)-x(1,:))));
ysp3=@(t) 2*mup(t).*mus(t).*beta(t).*t.*  exp(1i*(-mup(t)*x(2,:)-mus(t)*y(2,:)+t*(y(1,:)-x(1,:))));
ysp4=@(t) 2*mup(t).*beta(t).*t.^2 .*      exp(1i*(-mup(t)*x(2,:)-mus(t)*y(2,:)+t*(y(1,:)-x(1,:))));

yps1=@(t) 2*mus(t).*beta(t).*t.^2 .*      exp(1i*(-mus(t)*x(2,:)-mup(t)*y(2,:)+t*(y(1,:)-x(1,:))));
yps2=@(t) -2*mup(t).*mus(t).*beta(t).*t.* exp(1i*(-mus(t)*x(2,:)-mup(t)*y(2,:)+t*(y(1,:)-x(1,:))));
yps3=@(t) -2*beta(t).*t.^3 .*             exp(1i*(-mus(t)*x(2,:)-mup(t)*y(2,:)+t*(y(1,:)-x(1,:))));
yps4=@(t) 2*mup(t).*beta(t).*t.^2 .*      exp(1i*(-mus(t)*x(2,:)-mup(t)*y(2,:)+t*(y(1,:)-x(1,:))));

Gc1=integral(@(t) forec(t).*(ys1(t)+yp1(t)+ysp1(t)+yps1(t)), -l, l, 'ArrayValued',true);
Gc2=integral(@(t) forec(t).*(ys2(t)+yp2(t)+ysp2(t)+yps2(t)), -l, l, 'ArrayValued',true);
Gc3=integral(@(t) forec(t).*(ys3(t)+yp3(t)+ysp3(t)+yps3(t)), -l, l, 'ArrayValued',true);
Gc4=integral(@(t) forec(t).*(ys4(t)+yp4(t)+ysp4(t)+yps4(t)), -l, l, 'ArrayValued',true);

G(1,:)=G1(1,:)-G2(1,:)+Gc1;
G(2,:)=G1(2,:)-G2(2,:)+Gc2;
G(3,:)=G1(3,:)-G2(3,:)+Gc3;
G(4,:)=G1(4,:)-G2(4,:)+Gc4;

