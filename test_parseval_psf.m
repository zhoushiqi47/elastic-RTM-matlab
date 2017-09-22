%% ---------------- Parameters Setting---------------------------------------------%%


tic
n_recv = 1201;
lamda=1/2;
mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

testpoint=[0;10];
receiver = zeros(2,n_recv);


R = 0;
a = 600;

receiver(1,:) = linspace(-a,a,n_recv);    receiver(2,:) =  R*ones(1,n_recv);
%theta_r = (0:n_recv-1)*2*pi/n_recv;
%receiver(1,:)   = R*cos(theta_r); receiver(2,:)   = R*sin(theta_r);
   
x=0;
z=10;
%load pointdata G;
%% point scarttering data%%
G=Elastic_GreenTensor_Thalf_SIP(omega,kp,ks,[testpoint(1)*ones(1,n_recv); testpoint(2)*ones(1,n_recv)],receiver);
Gs=Elastic_GreenTensor_Thalf_SIP_ks(omega,kp,ks,[testpoint(1)*ones(1,n_recv); testpoint(2)*ones(1,n_recv)],receiver);
Gp=Elastic_GreenTensor_Thalf_SIP_kp(omega,kp,ks,[testpoint(1)*ones(1,n_recv); testpoint(2)*ones(1,n_recv)],receiver);

%% 回传点源数据 %%
 
          
            %gr=Traction_Dhalf(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);

           gr=TractionGreenTensor_2D(omega,kp,ks,[x*ones(1,n_recv); z*ones(1,n_recv)],receiver);
           grs=TractionGreenTensor_2D_ks(omega,kp,ks,[x*ones(1,n_recv); z*ones(1,n_recv)],receiver);
           grp=TractionGreenTensor_2D_kp(omega,kp,ks,[x*ones(1,n_recv); z*ones(1,n_recv)],receiver);
           
           %gr=Elastic_GreenTensor_2D(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
           
           G11 =trapz( receiver(1,:),(gr(1,:).*conj(G(1,:))+ gr(2,:).*conj(G(2,:))));
           Gss11 =trapz(  receiver(1,:),(grs(1,:).*conj(Gs(1,:))+ grs(2,:).*conj(Gs(2,:))));
           Gpp11 =trapz(  receiver(1,:),(grp(1,:).*conj(Gp(1,:))+ grp(2,:).*conj(Gp(2,:))));
           Gps11 =trapz(  receiver(1,:),(grp(1,:).*conj(Gs(1,:))+ grp(2,:).*conj(Gs(2,:))));
           Gsp11 =trapz(  receiver(1,:),(grs(1,:).*conj(Gp(1,:))+ grs(2,:).*conj(Gp(2,:))));
           G011=Gss11+Gpp11+Gps11+Gsp11;


        

toc