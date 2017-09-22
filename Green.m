function val=Green(k,z,z0)
z1=z(1,:)-z0(1,:);
z2=z(2,:)-z0(2,:);
d=sqrt(z1.*z1+z2.*z2);
val = 1i/4*besselh(0,k*d);