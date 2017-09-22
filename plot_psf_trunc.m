clear;


lamda=1;
mu=1;
omega= pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);
tt=linspace(0,pi,101);
t=-kp*cos(tt);
y=[0;10];
z=[-1;10];
mus=(ks^2-t.^2).^(1/2);
mup=(kp^2-t.^2).^(1/2);
k1=kp*(1+0.0000001i);
k2=ks*(1+0.0000001i);


beta=ks^2-2*t.^2;
gamma=t.^2+mus.*mup;
delta=beta.^2+4*t.^2.*mus.*mup;
Eys=exp(-1i*mus*(y(2))+1i*t*y(1));
Ezp=exp( 1i*mup*(z(2))-1i*t*z(1));
%coco=cos(mup*z(2)-mus*y(2)+t*(y(1)-z(1)));
fore=-1i.*mup./(mu*gamma.*delta*2*pi);
EE=Eys.*Ezp;
B=t.*(ks*ks-2*gamma).*t.*mus;
B1=t.*(ks*ks).*t.*mus.*(Eys.*Ezp);
B2=t.*(-2*gamma).*t.*mus.*(Eys.*Ezp);
f1=fore.*B.*EE;
f2=(-1i*mus*(y(2))+1i*t*y(1)+1i*mup*(z(2))-1i*t*z(1))/ks;
f3=1i./delta;
figure;
plot(tt,imag(f1))
figure;
plot(tt,imag(f2))
%figure;
%plot(tt,imag(f3))


G=integral(@(t) A3(t,z(:,:),y(:,:),k1,k2,mu) ,-kp,kp, 'ArrayValued',true)


