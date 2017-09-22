%% ---------------- Parameters Setting---------------------------------------------%%￥￥
%%2016 12 21 by shiqi
tic
n_recv = 401;
lamda=1;
mu=1;
omega= 2*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

testpoint=[0;10];
receiver = zeros(2,n_recv);


R = 0;
a = 200;

receiver(1,:) = linspace(-a,a,n_recv);    receiver(2,:) =  R*ones(1,n_recv);
%theta_r = (0:n_recv-1)*2*pi/n_recv;
%receiver(1,:)   = R*cos(theta_r); receiver(2,:)   = R*sin(theta_r);
   
Nx = 201;
Nz = 201;
NN = Nx*Nz;
x = linspace( -2, 2, Nx);
z = linspace( 8,12 , Nz);

%% point scarttering data%%
G=Elastic_GreenTensor_Thalf_SIP2(omega,kp,ks,[testpoint(1)*ones(1,n_recv); testpoint(2)*ones(1,n_recv)],receiver);
load TractionData_401_201 TractionData
%%把样本点展开成向量，略去两次循环%%
re = [repmat(x,1,Nz);reshape(repmat(z,Nx,1),1,NN)];
sample = reshape(repmat(re,n_recv,1),2,n_recv*NN);
G11=zeros(Nx,Nz);
G21=zeros(Nx,Nz);
G12=zeros(Nx,Nz);
G22=zeros(Nx,Nz);
%% 回传点源数据 %%
scaling = (receiver(1,end)-receiver(1,1))/n_recv;
 for jj=1:Nz         
            
           %Gr=TractionGreenTensor_2D(omega,kp,ks,sample,repmat(receiver,1,NN));
           %Gr=TractionDGreenTensor_2D2(omega,kp,ks,sample,repmat(receiver,1,NN));
           Gr=TractionData(:,:,jj);
           gr = zeros(Nx,n_recv,4);
           for j=1:4 
               gr(:,:,j) = reshape(Gr(j,:),n_recv,Nx).';
           end
           
           %gr=Elastic_GreenTensor_2D(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
           V1 = gr(:,:,1)*conj(G(1,:).')+gr(:,:,2)*conj(G(2,:).');
           V2 = gr(:,:,3)*conj(G(1,:).')+gr(:,:,4)*conj(G(2,:).');
           V3 = gr(:,:,1)*conj(G(3,:).')+gr(:,:,2)*conj(G(4,:).');
           V4 = gr(:,:,3)*conj(G(3,:).')+gr(:,:,4)*conj(G(4,:).');
           G11(:,jj)=scaling*V1;
           G21(:,jj)=scaling*V2;
           G12(:,jj)=scaling*V3;
           G22(:,jj)=scaling*V4;
 end

%% plot %%
[xx,zz]=meshgrid(x,z);
%save pointsf G11 G22 G12 G21

%figure, imagesc(x,z,real(G11)');colorbar;
figure
subplot(2,2,1)
 imagesc(x,z,-imag(G11)');colorbar;
colormap(jet)

subplot(2,2,2)
 imagesc(x,z,-imag(G12)');colorbar;
colormap(jet)

subplot(2,2,3)
 imagesc(x,z,-imag(G21)');colorbar;
colormap(jet)

subplot(2,2,4)
 imagesc(x,z,-imag(G22)');colorbar;
colormap(jet)

figure
subplot(2,2,1)
 imagesc(x,z,real(G11)');colorbar;
colormap(jet)

subplot(2,2,2)
 imagesc(x,z,real(G12)');colorbar;
colormap(jet)

subplot(2,2,3)
 imagesc(x,z,real(G21)');colorbar;
colormap(jet)

subplot(2,2,4)
 imagesc(x,z,real(G22)');colorbar;
colormap(jet)

%figure, mesh(xx,zz,real(G11)');
%figure, mesh(xx,zz,-imag(G11)');
toc