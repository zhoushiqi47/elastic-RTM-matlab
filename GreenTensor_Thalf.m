function G=GreenTensor_Thalf(omega,kp,ks,x,y)
%x is point source
w2 = 1/omega^2 ;
k1=kp*(1+0.00000001i);
k2=ks*(1+0.00000001i);

n=size(x);
G=zeros(4,n(2));
l=5*(kp+ks);
xx=x;
xx(2,:)=-x(2,:);
G1=Elastic_GreenTensor_2D(omega,kp,ks,x,y);
G2=Elastic_GreenTensor_2D(omega,kp,ks,xx,y);
ndots=10000;
t=linspace(-l,l,ndots+1).';

d=2*l/ndots;

mup=(k1*k1-t.*t).^(1/2);
mus=(k2*k2-t.*t).^(1/2);
beta=k2*k2-2*t.*t;
delta=beta.*beta+4*t.*t.*mus.*mup;
forec=1i*w2./(2*pi*repmat(delta,1,n(2)));

ys1=repmat(mus.*beta.*beta,1,n(2)).*      exp(1i*(-mus*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
ys2=repmat(t.*beta.*beta,1,n(2)).*        exp(1i*(-mus*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
ys3=repmat(4*mus.*mup.*t.^3,1,n(2)).*     exp(1i*(-mus*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
ys4=repmat(4*mup.*t.^4,1,n(2)).*          exp(1i*(-mus*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));

yp1=repmat(4*mus.*t.^4,1,n(2)).*          exp(1i*(-mup*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
yp2=repmat(-4*mus.*mup.*t.^3,1,n(2)).*    exp(1i*(-mup*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
yp3=repmat(-t.*beta.*beta,1,n(2)).*       exp(1i*(-mup*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));
yp4=repmat(mup.*beta.*beta,1,n(2)).*      exp(1i*(-mup*(x(2,:)+y(2,:))+t*(y(1,:)-x(1,:))));

ysp1=repmat(2*mus.*beta.*t.^2,1,n(2)).*   exp(1i*(-mus*x(2,:)-mup*y(2,:)+t*(y(1,:)-x(1,:))));
ysp2=repmat(2*beta.*t.^3,1,n(2)).*        exp(1i*(-mus*x(2,:)-mup*y(2,:)+t*(y(1,:)-x(1,:))));
ysp3=repmat(2*mup.*mus.*beta.*t,1,n(2)).* exp(1i*(-mus*x(2,:)-mup*y(2,:)+t*(y(1,:)-x(1,:))));
ysp4=repmat(2*beta.*mup.*t.^2,1,n(2)).*   exp(1i*(-mus*x(2,:)-mup*y(2,:)+t*(y(1,:)-x(1,:))));

yps1=repmat(2*beta.*mus.*t.^2,1,n(2)).*   exp(1i*(-mup*x(2,:)-mus*y(2,:)+t*(y(1,:)-x(1,:))));
yps2=repmat(-2*mup.*mus.*beta.*t,1,n(2)).*exp(1i*(-mup*x(2,:)-mus*y(2,:)+t*(y(1,:)-x(1,:))));
yps3=repmat(-2*beta.*t.^3,1,n(2)).*       exp(1i*(-mup*x(2,:)-mus*y(2,:)+t*(y(1,:)-x(1,:))));
yps4=repmat(2*mup.*beta.*t.^2,1,n(2)).*   exp(1i*(-mup*x(2,:)-mus*y(2,:)+t*(y(1,:)-x(1,:))));

Gc1=trapz(forec.*(ys1+yp1+ysp1+yps1),t);
Gc2=trapz(forec.*(ys2+yp2+ysp2+yps2),t);
Gc3=trapz(forec.*(ys3+yp3+ysp3+yps3),t);
Gc4=trapz(forec.*(ys4+yp4+ysp4+yps4),t);

G(1,:)=G1(1,:)-G2(1,:)+Gc1;
G(2,:)=G1(2,:)-G2(2,:)+Gc2;
G(3,:)=G1(3,:)-G2(3,:)+Gc3;
G(4,:)=G1(4,:)-G2(4,:)+Gc4;