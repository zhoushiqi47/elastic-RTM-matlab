function f=A3(t,z,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
gamma=t^2+mus*mup;
delta=beta^2+4*t^2*mus*mup;
Eys=exp(-1i*mus*(y(2,:))+1i*t*y(1,:));
Eyp=exp(-1i*mup*(y(2,:))+1i*t*y(1,:));
Ezs=exp( 1i*mus*(z(2,:))-1i*t*z(1,:));
Ezp=exp( 1i*mup*(z(2,:))-1i*t*z(1,:));
fore=-1i/(mu*gamma*delta*2*pi);
%B=t*(ks*ks-2*gamma)*t*mus*(Eys.*Ezp);
B=exp(-1i*mus*(y(2,:))+1i*mup*(z(2,:)));
f=fore*(B);
end