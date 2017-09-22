function G=TractionGreenTensor_2D_kp(omega,kp,ks,z0,z)
%% n=(0,1)
%% operator T to z
%% compute the traction green tensor of elastic wave operator in 2D
%% shear wave velocity

%% pressure wave Green Tensor
w2 = 1/omega^2 ;
mu = omega.^2./(ks.^2);
cp2 = omega.^2./(kp.^2);
lamda = cp2 - 2*mu;
%% the second method for implementing the elastic green's tensor

z1=z(1,:)-z0(1,:);
z2=z(2,:)-z0(2,:);
r=sqrt(z1.*z1+z2.*z2);

Hp0 = besselh(0,kp*r);
Hp1 = besselh(1,kp*r);


r2 = r.*r;
r3= r2.*r;


phip1 = 1i/4*w2*(-kp^2*Hp0 + 2*kp*Hp1./r);


phip2 = 1i/4*w2*(kp^3*Hp1 + 4*kp^2*Hp0./r - 8*kp*Hp1./r2);



gp111 = 3*phip1.*z1./r2 + phip2.*z1.^3./r3;
gp112 = phip1.*z2./r2 + phip2.*z1.^2.*z2./r3;
gp122 = phip1.*z1./r2 + phip2.*z2.^2.*z1./r3;
gp222 = 3*phip1.*z2./r2 + phip2.*z2.^3./r3;

G(1,:) =  - 2*mu*gp112;
G(2,:) =  - cp2*gp122 - lamda*gp111;
G(3,:) =  - 2*mu*gp122;
G(4,:) =  - cp2*gp222 - lamda*gp112;



