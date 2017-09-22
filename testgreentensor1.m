tic
n_recv = 201;
lamda=1/2;
mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

testpoint=[0;-10];

Nx = 51;
Nz = 51;
N = Nx*Nz;
x = linspace( -2, 2, Nx);
z = linspace( -12,-8 , Nz);
G=zeros(Nx,Nz,4);

for j=1:Nx
    for k=1:Nz
        G(j,k,:)=Elastic_GreenTensor_Thalf(omega,kp,ks,testpoint,[x(j);z(k)]);
    end
end


          
figure
subplot(2,2,1)
 imagesc(x,z,-imag(G(:,:,1))');colorbar;
colormap(jet)

subplot(2,2,2)
 imagesc(x,z,-imag(G(:,:,2))');colorbar;
colormap(jet)

subplot(2,2,3)
 imagesc(x,z,-imag(G(:,:,3))');colorbar;
colormap(jet)

subplot(2,2,4)
 imagesc(x,z,-imag(G(:,:,4))');colorbar;
colormap(jet)           

%figure, mesh(xx,zz,real(G(:,:,1))');
toc