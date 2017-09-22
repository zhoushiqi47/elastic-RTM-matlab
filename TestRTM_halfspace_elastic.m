%% ---------------- Parameters Setting---------------------------------------------%%

clear;
n_src = 401;
n_recv = 401;

lamda=1;
mu=1;
omega= 2*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

npts = 256 ; %% discretization point of integral equation

source   = zeros(2,n_src);
receiver = zeros(2,n_recv);

Q=[1 ,0;0 ,1];
R = 0;
a = 200;

%% haf space reflection
source(1,:)   = linspace(-a,a,n_src);     source(2,:)   =  R*ones(1,n_src);
receiver(1,:) = linspace(-a,a,n_recv);    receiver(2,:) =  R*ones(1,n_recv);
bctype = 1;  %% 1 circle;2 kite;            


tic

 %% the obstacle is under the half space
 % 合成散射数据
 %Data = NyElasticWave_src_half(omega , lamda, mu, npts, bctype, n_src, n_recv, source, receiver);
 %save scatteredfield_256_401_r1_h10 Data   
 load scatteredfield_256_401_r1_h10 Data  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Add additive Gaussian random noise.
%     data(:,:,k) = U(:,:,k);
%     maxU = max(max(abs( U(:,:,k))));
%     noiselevel = 0.2*maxU;
 %  U(:,:,k) =  U(:,:,k) + noiselevel/sqrt(2)*maxU*((2*rand(size( U(:,:,k)))-1)+i*(2*rand(size( U(:,:,k)))-1));






%% Sampling domain;
Nx = 201;
Nz = 201;
x = linspace( -2, 2, Nx);
z = linspace( 8, 12, Nz);

%% 逆时偏移成像
Irtm = RTM_halfspace_elastic2(n_src, n_recv,source, receiver, x, z , Data, omega , kp ,ks);
%load circle_r1_401_-10_256_201 Irtm


[xx,zz]=meshgrid(x,z);
 
figure, imagesc(x,z,imag(Irtm)');colorbar; 
colormap(jet)

%% Plot the true scatter
n = npts;
node = 0:2*n-1;
t = pi*node(:)/n;
[xt,zt]=circlebc(t,1);
hold on;
plot(xt,zt,'r');

figure, imagesc(x,z,real(Irtm)');colorbar; 
colormap(jet)

hold on;
plot(xt,zt,'r');

