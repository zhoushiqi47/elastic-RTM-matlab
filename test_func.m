

clear
clc

nx = 201;
nz = 201;
x=linspace(-4,4,nx);
z=linspace(100,104,nz);
x0=1;
z0=103;

k=2*pi*1;
N=127;
theta = pi/N*(0:N-1);
F = zeros(nx, nz);
F2 = zeros(nx, nz);
n_recv = 128;
theta_r = (0:n_recv-1)*pi/(n_recv);


flag=1;  %% flag=1 near fieds flag=2, far field
if flag==1
    rho = 100000;
else 
    rho=1;
end
receiver(1,:) = rho*cos(theta_r); receiver(2,:) = rho*sin(theta_r);
a = 20; rho=8 ;
 
receiver(1,:) = linspace(-a,a,n_recv);receiver(2,:) = 0*ones(1,n_recv);

g2 =  Green(k,[x0;z0]*ones(1,n_recv),receiver);
Gn2 = Green_Grad(k, receiver, [x0;z0]*ones(1,n_recv));
Gn4 = Green_Grad(k, receiver, [x0;-z0]*ones(1,n_recv));
g4 =  Green(k,[x0;-z0]*ones(1,n_recv),receiver);
tic
for ix=1:nx
    for iz=1:nz
       
        t1 =  x(ix)*cos(theta) + z(iz)*sin(theta);
        t2 =  x(ix)*cos(theta) - z(iz)*sin(theta);
        t3 =  x0*cos(theta) + z0*sin(theta);
        t4 =  x0*cos(theta) - z0*sin(theta);
        F(ix,iz) = pi/N*sum( ( exp(sqrt(-1)*k*t1)-exp(sqrt(-1)*k*t2)).*( exp(-sqrt(-1)*k*t3)-exp(-sqrt(-1)*k*t4)));
        g1 = Green(k,[x(ix);z(iz)]*ones(1,n_recv),receiver);
        g3 = Green(k,[x(ix);-z(iz)]*ones(1,n_recv),receiver);
        F2(ix,iz) = k*rho*pi/n_recv*sum(( Gn2(2,:)   ).*( conj( g1 )  ) );
    end
end
toc
figure
 [xx,zz]=meshgrid(x,z);
subplot(1,2,1);
mesh(xx,zz,real(F)');
subplot(1,2,2);
mesh(xx,zz,imag(F)');
[xx,zz]=meshgrid(x,z);
p = 1/4*besselj(0,k*sqrt(abs(xx-x0).^2 +abs(z0-zz).^2));

figure

subplot(1,2,1);
mesh(xx,zz,real(F2)');
subplot(1,2,2);
mesh(xx,zz,imag(F2)');

