function Irtm=RTM_halfspace_elastic2(n_src, n_recv,source, receiver, x, z , Data, omegas , kp ,ks)
%% Reverse Time Migration for obstacle scattering data available on a line.
%% n_src, n_recv denotes the number of source and receiver
%% Data(n_recv,n_src,kp，ks)
%% 2016 7 23 shiqi
tic
N = length(omegas);
Nx = length(x);
Nz = length(z);


    %% 半空间逆时偏移成像
scaling = 1/N*(receiver(1,end)-receiver(1,1))*(source(1,end)-source(1,1))/n_src/n_recv;

%%把样本点展开成向量，略去两次循环
%%当Nx Nz 比较大的时候 有可能会内纯溢出 此时只用一次循环



irtm = zeros(Nx,Nz);
Irtm = zeros(Nx,Nz);
load TractionData_401_201 TractionData
for k=1:N
    omega = omegas(k);
    
    for jj=1:Nz
        gr = zeros(Nx,n_recv,4);
        gs = zeros(Nx,n_src,4);
        re=[x;repmat(z(jj),1,Nx)];
        sample1 = reshape(repmat(re,n_recv,1),2,n_recv*Nx);
        sample2 = reshape(repmat(re,n_src,1),2,n_src*Nx);
        
    if n_src==n_recv
        %Gr=TractionGreenTensor_2D(omega,kp,ks,sample1,repmat(receiver,1,Nx));
        Gr=TractionData(:,:,jj);
        Gs=Gr;
    else
        Gr=TractionDGreenTensor_2D2(omega,kp,ks,sample1,repmat(receiver,1,Nx));
        Gs=TractionDGreenTensor_2D2(omega,kp,ks,sample2,repmat(source,1,Nx));
    end
    
    for j=1:4
        gr(:,:,j) = reshape(Gr(j,:),n_recv,Nx).';
        gs(:,:,j) = reshape(Gs(j,:),n_src,Nx).';
    end
    %水平和垂直两个方向激发后叠加
    %% 反传共轭散射数据
    Vs1 = gr(:,:,1)*conj(Data(:,:,1,1,k))+gr(:,:,2)*conj(Data(:,:,2,1,k));
    Vs2 = gr(:,:,3)*conj(Data(:,:,1,1,k))+gr(:,:,4)*conj(Data(:,:,2,1,k));
    Vs3 = gr(:,:,1)*conj(Data(:,:,1,2,k))+gr(:,:,2)*conj(Data(:,:,2,2,k));
    Vs4 = gr(:,:,3)*conj(Data(:,:,1,2,k))+gr(:,:,4)*conj(Data(:,:,2,2,k));
    % 和相对应的点源做互相关
    irtm(:,jj) =  scaling*sum(Vs1.*gs(:,:,1) +Vs2.*gs(:,:,3) + Vs3.*gs(:,:,2) + Vs4.*gs(:,:,4),2);
    end

    Irtm  = Irtm + irtm;
end

disp('半空间逆时偏移结束');
toc
end



    
            