function G=Elastic_GreenTensor_2D(omega,kp,ks,z,z0)

%% compute the dynaic green tensor of elastic wave operator in 2D
%% shear wave velocity
cs = omega/ks;
%% pressure wave Green Tensor
w2 = 1/omega^2 ;
% Gp = -w2*Green_DivGrad_2D(kp,z,z0) ;
% 
% %% shear wave Green Tensor;
% gs = Green(ks,z,z0);
% Gs = w2*Green_DivGrad_2D(ks,z,z0) ;
% Gs(1,:) = Gs(1,:) + 1/cs^2 * gs;
% Gs(4,:) = Gs(4,:) + 1/cs^2 * gs;
% 
% %% the whole Green's Tensor G = Gp + Gs;
% 
% G = Gp + Gs;


%% the second method for implementing the elastic green's tensor

z1=z(1,:)-z0(1,:);
z2=z(2,:)-z0(2,:);
r=sqrt(z1.*z1+z2.*z2);

Hp0 = besselh(0,kp*r);
Hp1 = besselh(1,kp*r);

Hs0 = besselh(0,ks*r);
Hs1 = besselh(1,ks*r);

r2 = r.*r;
phi1 = 1i/4*( 1/cs^2*Hs0 - w2*( (ks*Hs1 - kp*Hp1)./r ) );
phi2 = 1i/4*w2*( 2*( (ks*Hs1 - kp*Hp1)./r ) - ks^2.*Hs0 + kp^2*Hp0) ;


G(1,:) = phi1 + phi2.*z1.*z1./r2;
G(2,:) =        phi2.*z1.*z2./r2;
G(3,:) = G(2,:);
G(4,:) = phi1 + phi2.*z2.*z2./r2;