tic

lamda=1/2;
mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

testpoint=[0;10];

Nx = 101;
Nz = 101;
N = Nx*Nz;
x = linspace( -2, 2, Nx);
z = linspace( 8,12 , Nz);
G=zeros(Nx,Nz,4);
%%把样本点展开成向量，略去两次循环
re=[repmat(x,1,Nz);reshape(repmat(z,Nx,1),1,N)];


          
           
reG=Elastic_GreenTensor_Thalf_SIP(omega,kp,ks,repmat(testpoint,1,N),re);
           
          
 
for j=1:4
    G(:,:,j)=reshape(reG(j,:),Nx,Nz);
end

figure
subplot(2,2,1)
 imagesc(x,z,imag(G(:,:,1))');colorbar;
colormap(jet)

subplot(2,2,2)
 imagesc(x,z,imag(G(:,:,2))');colorbar;
colormap(jet)

subplot(2,2,3)
 imagesc(x,z,imag(G(:,:,3))');colorbar;
colormap(jet)

subplot(2,2,4)
 imagesc(x,z,imag(G(:,:,4))');colorbar;
colormap(jet)  
%figure, mesh(xx,zz,real(G(:,:,1))');
toc