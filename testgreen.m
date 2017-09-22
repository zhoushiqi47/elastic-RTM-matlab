tic

lamda=1/2;
mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

testpoint=[0;-10];

Nx = 101;
Nz = 101;
N = Nx*Nz;
x = linspace( -2, 2, Nx);
z = linspace( -12,-8 , Nz);

%%把样本点展开成向量，略去两次循环
re=[repmat(x,1,Nz);reshape(repmat(z,Nx,1),1,N)];


          
           
G=Green(omega,repmat(testpoint,1,N),re);
 reG=reshape(G,Nx,Nz);          
          


figure

 imagesc(x,z,imag(reG.'));colorbar;

colormap(jet)  
%figure, mesh(xx,zz,real(G(:,:,1))');
toc