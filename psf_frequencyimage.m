%% ---------------- Parameters Setting---------------------------------------------%%￥￥
%%2016 12 22 by shiqi
lamda=1;
mu=1;
omega= 2*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

testpoint=[0;10];
Nx = 201;
Nz = 201;
NN = Nx*Nz;
x = linspace( -2, 2, Nx);
z = linspace( 8, 12 , Nz);
re = [repmat(x,1,Nz);reshape(repmat(z,Nx,1),1,NN)];
%% 计算样本点上psf的值 是一个长向量
[psflong1,psflong2]=psf_frequency(omega,kp,ks,re,repmat(testpoint,1,NN));

%% 变回矩阵
psf1= reshape(psflong1,Nx,Nz);

psf2= reshape(psflong2,Nx,Nz);
figure
subplot(2,2,1)
 imagesc(x,z,-imag(psf1)');colorbar;
colormap(jet)

subplot(2,2,2)
 imagesc(x,z,-imag(psf2)');colorbar;
colormap(jet)

subplot(2,2,3)
 imagesc(x,z,real(psf1)');colorbar;
colormap(jet)

subplot(2,2,4)
 imagesc(x,z,real(psf2)');colorbar;
colormap(jet)

%figure
%imagesc(x,z,-imag(psf1+psf2)');colorbar
