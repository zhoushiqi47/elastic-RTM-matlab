function Irtm=RTM_halfspace_elastic(n_src, n_recv,source, receiver, x, z , Data, omegas , kp ,ks)
%% Reverse Time Migration for obstacle scattering data available on a line.
%% n_src, n_recv denotes the number of source and receiver
%% Data(n_recv,n_src,kp，ks)
%% 2016 7 23 shiqi
tic
N = length(omegas);
Nx = length(x);
Nz = length(z);
NN = Nx*Nz;


    %% 半空间逆时偏移成像
scaling = 1/N*(receiver(1,end)-receiver(1,1))*(source(1,end)-source(1,1))/n_src/n_recv;

%%把样本点展开成向量，略去两次循环
re = [repmat(x,1,Nz);reshape(repmat(z,Nx,1),1,NN)];
sample1 = reshape(repmat(re,n_recv,1),2,n_recv*NN);
sample2 = reshape(repmat(re,n_src,1),2,n_recv*NN);
gr = zeros(NN,n_recv,4);
gs = zeros(NN,n_src,4);
p = zeros(NN,1);
for k=1:N
    omega = omegas(k);
    if n_src==n_recv
        Gr=TractionDGreenTensor_2D2(omega,kp,ks,sample1,repmat(receiver,1,NN));
        Gs=Gr;
    else
        Gr=TractionDGreenTensor_2D2(omega,kp,ks,sample1,repmat(receiver,1,NN));
        Gs=TractionDGreenTensor_2D2(omega,kp,ks,sample2,repmat(source,1,NN));
    end
    
    for j=1:4
        gr(:,:,j) = reshape(Gr(j,:),n_recv,NN).';
        gs(:,:,j) = reshape(Gs(j,:),n_src,NN).';
    end
    save Dincidentfield_201_10_2 gr gs
    %load incidentfield_201_-10_1 gr gs
    %水平和垂直两个方向激发后叠加
    %% 反传共轭散射数据
    Vs1 = gr(:,:,1)*conj(Data(:,:,1,1,k))+gr(:,:,2)*conj(Data(:,:,2,1,k));
    Vs2 = gr(:,:,3)*conj(Data(:,:,1,1,k))+gr(:,:,4)*conj(Data(:,:,2,1,k));
    Vs3 = gr(:,:,1)*conj(Data(:,:,1,2,k))+gr(:,:,2)*conj(Data(:,:,2,2,k));
    Vs4 = gr(:,:,3)*conj(Data(:,:,1,2,k))+gr(:,:,4)*conj(Data(:,:,2,2,k));
    % 和相对应的点源做互相关
    p = p + scaling*sum(Vs1.*gs(:,:,1) +Vs2.*gs(:,:,3) + Vs3.*gs(:,:,2) + Vs4.*gs(:,:,4),2);
end

Irtm  = reshape(p,Nx,Nz);

disp('半空间逆时偏移结束');
toc
end



    
            