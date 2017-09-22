%%提前计算好 反传的traction数据 并以z方向排列，所有x方向在同一个向量


n_src = 401;
n_recv = 401;
R = 0;
a = 200;
receiver = zeros(2,n_recv);
receiver(1,:) = linspace(-a,a,n_recv);    receiver(2,:) =  R*ones(1,n_recv);

lamda=1;
mu=1;
omega= 2*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

Nx = 201;
Nz = 201;
x = linspace( -2, 2, Nx);
z = linspace( 8, 12, Nz);

TractionData=zeros(4,n_recv*Nx,Nz);

 parfor jj=1:Nz
        re=[x;repmat(z(jj),1,Nx)];
        sample1 = reshape(repmat(re,n_recv,1),2,n_recv*Nx);
        Gr=TractionDGreenTensor_2D2(omega,kp,ks,sample1,repmat(receiver,1,Nx));
        TractionData(:,:,jj)=Gr;
 end
        
        