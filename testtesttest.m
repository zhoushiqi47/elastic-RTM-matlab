%% ---------------- Parameters Setting---------------------------------------------%%

clear;
n_src = 3;
n_recv =3;

lamda=1/2;
mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

npts = 64 ; %% discretization point of integral equation

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
 Data = ElasticWave_src_half(omega , lamda, mu, npts, bctype, n_src, n_recv, source, receiver, Q);