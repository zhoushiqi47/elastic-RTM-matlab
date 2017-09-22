function Irtm=RTM_halfspace_elastic1(n_src, n_recv,source, receiver, x, z , Data, omega , kp ,ks)
%% Reverse Time Migration for obstacle scattering data available on a line.
%% n_src, n_recv denotes the number of source and receiver
%% Data(n_recv,n_src,kp，ks)

N = length(omega);
Nx = length(x);
Nz = length(z);

Irtm = zeros(Nx,Nz);
tic

    %% 半空间逆时偏移成像
scaling = 1/N*(receiver(1,end)-receiver(1,1))*(source(1,end)-source(1,1))/n_src/n_recv;

for ix=  1:Nx
    for iz=1:Nz
        for k=1:N

            if n_src==n_recv
               Gr=TractionDGreenTensor_2D2(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
               Gs=Gr;
            else
               Gr=TractionDGreenTensor_2D2(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
               Gs=TractionDGreenTensor_2D2(omega,kp,ks,[x(ix)*ones(1,n_src); z(iz)*ones(1,n_src)],source);
            end
            
            %水平和垂直两个方向激发后叠加
           %% 反传共轭散射数据
            Vs1 = Gr(1,:)*conj(Data(:,:,1,1,k))+Gr(2,:)*conj(Data(:,:,2,1,k));
            Vs2 = Gr(3,:)*conj(Data(:,:,1,1,k))+Gr(4,:)*conj(Data(:,:,2,1,k));
            Vs3 = Gr(1,:)*conj(Data(:,:,1,2,k))+Gr(2,:)*conj(Data(:,:,2,2,k));
            Vs4 = Gr(3,:)*conj(Data(:,:,1,2,k))+Gr(4,:)*conj(Data(:,:,2,2,k));
            
            % 和相对应的点源做互相关

            Irtm(ix,iz)=Irtm(ix,iz) + scaling*sum(Gs(1,:).*Vs1 + Gs(3,:).*Vs2 + Gs(2,:).*Vs3 + Gs(4,:).*Vs4);
      

        end
    end
end

disp('半空间逆时偏移结束');
toc
end



    
            