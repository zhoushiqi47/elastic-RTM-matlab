function [G1,G2]=psf_frequency(omega,kp,ks,z,y)
%% computate psf in frequency field by paserval theorem
%% y is the point soure,z is the sample point
%% the integral path is only on (-kp,kp)

% 向量形式积分
%y is point source
% 2016 12 22 by shiqi



mu = omega.^2./(ks.^2);
l=3*(kp+ks);
k1=kp*(1+0.0000001i);
k2=ks*(1+0.0000001i);

G11=integral(@(t) B21(t,z(:,:),y(:,:),k1,k2,mu) ,-kp,kp, 'ArrayValued',true);
G12=integral(@(t) B22(t,z(:,:),y(:,:),k1,k2,mu) ,-ks,ks, 'ArrayValued',true);
G1=G11+G12;
%G1=integral(@(t) B1(t,z(:,:),y(:,:),k1,k2,mu) ,-l,l, 'ArrayValued',true);

G2=integral(@(t) B2(t,z(:,:),y(:,:),k1,k2,mu) ,-ks,ks, 'ArrayValued',true);
end

%%integral on【-kp,kp】

function f=A1(t,z,y,kp,ks,mu)
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
A=ks*ks*mus*(mup*mus*Eys.*Ezs+t*t*Eyp.*Ezp);
B=t*(ks*ks-2*gamma)*t*mus*(Eys.*Ezp-Eyp.*Ezs);
f=fore*(A+B);
end

function f=A2(t,z,y,kp,ks,mu)
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
A=ks*ks*mus*(mup*mus*Eys.*Ezs+t*t*Eyp.*Ezp);
f=fore*(A);
end

function f=A22(t,z,y,kp,ks,mu)
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
A=ks*ks*mus*mup*mus*Eys.*Ezs;
f=fore*(A);
end

function f=A21(t,z,y,kp,ks,mu)
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
A=ks*ks*mus*(t*t*Eyp.*Ezp);
f=fore*(A);
end

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
B=t*(ks*ks-2*gamma)*t*mus*(Eys.*Ezp);
f=fore*(B);
end
function f=A31(t,z,y,kp,ks,mu)
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
B=t*(ks*ks)*t*mus*(Eys.*Ezp);
f=fore*(B);
end
function f=A32(t,z,y,kp,ks,mu)
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
B=t*(-2*gamma)*t*mus*(Eys.*Ezp);
f=fore*(B);
end

function f=A4(t,z,y,kp,ks,mu)
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
B=t*(ks*ks-2*gamma)*t*mus*(-Eyp.*Ezs);
f=fore*(B);
end
%% integral on [-infty, infty]
function f=B1(t,z,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
gamma=t^2+mus*mup;
delta=beta^2+4*t^2*mus*mup;
Ts=[mus*mup, t*mus];
Tp=[t*t, -t*mus];
Ns=conj([mus*beta; 2*t*mus*mup]);
Np=conj([2*t*t*mus; -t*beta]);
Eys=conj(exp( 1i*mus*(y(2,:))-1i*t*y(1,:)));
Eyp=conj(exp( 1i*mup*(y(2,:))-1i*t*y(1,:)));
Ezs=exp( 1i*mus*(z(2,:))-1i*t*z(1,:));
Ezp=exp( 1i*mup*(z(2,:))-1i*t*z(1,:));
fore=-1i/(mu*gamma*conj(delta)*2*pi);
A=Ts*Ns*Eys.*Ezs+Tp*Np*Eyp.*Ezp;
B=Tp*Ns*Eys.*Ezp+Ts*Np*Eyp.*Ezs;
f=fore*(A+B);
end

function f=B2(t,z,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
gamma=t^2+mus*mup;
delta=beta^2+4*t^2*mus*mup;
Ts=[mus*mup, t*mus];
Tp=[t*t, -t*mus];
Ns=conj([mus*beta; 2*t*mus*mup]);
Np=conj([2*t*t*mus; -t*beta]);
Eys=conj(exp( 1i*mus*(y(2,:))-1i*t*y(1,:)));
Eyp=conj(exp( 1i*mup*(y(2,:))-1i*t*y(1,:)));
Ezs=exp( 1i*mus*(z(2,:))-1i*t*z(1,:));
Ezp=exp( 1i*mup*(z(2,:))-1i*t*z(1,:));
fore=-1i/(mu*gamma*conj(delta)*2*pi);
A=Ts*Ns*Eys.*Ezs+Tp*Np*Eyp.*Ezp;
f=fore*(A);
end

function f=B21(t,z,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
gamma=t^2+mus*mup;
delta=beta^2+4*t^2*mus*mup;
Ts=[mus*mup, t*mus];
Tp=[t*t, -t*mus];
Ns=conj([mus*beta; 2*t*mus*mup]);
Np=conj([2*t*t*mus; -t*beta]);
Eys=conj(exp( 1i*mus*(y(2,:))-1i*t*y(1,:)));
Eyp=conj(exp( 1i*mup*(y(2,:))-1i*t*y(1,:)));
Ezs=exp( 1i*mus*(z(2,:))-1i*t*z(1,:));
Ezp=exp( 1i*mup*(z(2,:))-1i*t*z(1,:));
fore=-1i/(mu*gamma*conj(delta)*2*pi);
%A=Ts*Ns*Eys.*Ezs;
A=Tp*Np*Eyp.*Ezp;
f=fore*(A);
end

function f=B22(t,z,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
gamma=t^2+mus*mup;
delta=beta^2+4*t^2*mus*mup;
Ts=[mus*mup, t*mus];
Tp=[t*t, -t*mus];
Ns=conj([mus*beta; 2*t*mus*mup]);
Np=conj([2*t*t*mus; -t*beta]);
Eys=conj(exp( 1i*mus*(y(2,:))-1i*t*y(1,:)));
Eyp=conj(exp( 1i*mup*(y(2,:))-1i*t*y(1,:)));
Ezs=exp( 1i*mus*(z(2,:))-1i*t*z(1,:));
Ezp=exp( 1i*mup*(z(2,:))-1i*t*z(1,:));
fore=-1i/(mu*gamma*conj(delta)*2*pi);
A=Ts*Ns*Eys.*Ezs;
%A=Tp*Np*Eyp.*Ezp;
f=fore*(A);
end

function f=B3(t,z,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
gamma=t^2+mus*mup;
delta=beta^2+4*t^2*mus*mup;
Ts=[mus*mup, t*mus];
Tp=[t*t, -t*mus];
Ns=conj([mus*beta; 2*t*mus*mup]);
Np=conj([2*t*t*mus; -t*beta]);
Eys=conj(exp( 1i*mus*(y(2,:))-1i*t*y(1,:)));
Eyp=conj(exp( 1i*mup*(y(2,:))-1i*t*y(1,:)));
Ezs=exp( 1i*mus*(z(2,:))-1i*t*z(1,:));
Ezp=exp( 1i*mup*(z(2,:))-1i*t*z(1,:));
fore=-1i/(mu*gamma*conj(delta)*2*pi);
B=Ts*Np*Eyp.*Ezs;
f=fore*(B);
end
