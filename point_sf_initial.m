%% ---------------- Parameters Setting---------------------------------------------%%


tic
n_recv = 401;
lamda=1/2;
mu=1/4;
omega=pi;
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
x = linspace( -2, 2, Nx);
z = linspace( 8,12 , Nz);
G11=zeros(Nx,Nz);
G21=zeros(Nx,Nz);
G12=zeros(Nx,Nz);
G22=zeros(Nx,Nz);
%load pointdata G;
%% point scarttering data%%
G=Elastic_GreenTensor_Thalf_SIP2(omega,kp,ks,[testpoint(1)*ones(1,n_recv); testpoint(2)*ones(1,n_recv)],receiver);

%% 回传点源数据 %%
  for ix=  1:Nx
        for iz=1:Nz
          
            %gr=Traction_Dhalf(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);

           gr=TractionGreenTensor_2D(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
           
           %gr=Elastic_GreenTensor_2D(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
           G11(ix,iz) =trapz(  receiver(1,:),(gr(1,:).*conj(G(1,:))+ gr(2,:).*conj(G(2,:))));
           G21(ix,iz) =trapz(  receiver(1,:),(gr(3,:).*conj(G(1,:))+ gr(4,:).*conj(G(2,:))));
           G12(ix,iz) =trapz(  receiver(1,:),(gr(1,:).*conj(G(3,:))+ gr(2,:).*conj(G(4,:))));
           G22(ix,iz) =trapz(  receiver(1,:),(gr(3,:).*conj(G(3,:))+ gr(4,:).*conj(G(4,:))));
        end
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

%figure, mesh(xx,zz,real(G11)');
%figure, mesh(xx,zz,-imag(G11)');
toc