clear;


n_src = 2 ;
n_recv = 11;

flag = 1;
theta_r = (0:n_recv-1)*2*pi/n_recv;
theta_s = (0:n_src-1)*2*pi/n_src ;

n = 64; %% discretization point of integral equation
pp = 1;

lamda=1/2;
mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

% for R=[10 8 6 4 3.5 3.2 3.1 2.99 2.49];
%source   = zeros(2,n_src);
receiver = zeros(2,n_recv);
R = 5;
a = 0;
% source(1,:)   = R*cos(theta_s); source(2,:)   = R*sin(theta_s);
source=[0 0;-10 -10.2];
%receiver(1,:) = R*cos(theta_r); receiver(2,:) = R*sin(theta_r)-10;
receiver(1,:) = linspace(-R,R,n_recv);    receiver(2,:) =  a*ones(1,n_recv);
Q=[1 ,0;0 ,1];
bctype = 1;
srctype=1;
%% solve via Nystrom's method

%% U =   NystromImpedance(npts, n_src, n_recv, bctype, freq, source, receiver);
Data = NyElasticWave_src_half1(omega, lamda, mu, n, bctype, n_src, n_recv, source, receiver);

 
Data2= Elastic_GreenTensor_Thalf1(omega,kp,ks,source(:,1)*ones(1,n_recv),receiver);

Data3= Elastic_GreenTensor_Thalf1(omega,kp,ks,source(:,2)*ones(1,n_recv),receiver);

[Data(:,:,1,1) -Data2(1,:).' -Data3(1,:).']
