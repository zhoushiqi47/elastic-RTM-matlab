function G=Traction_Dhalf1(omega,kp,ks,x,y)
%% traction tensor of displacement free green tensor
%缺少表面波部分
%x is point source
xx=x;
xx(2,:)=-x(2,:);
G1=TractionGreenTensor_2D(omega,kp,ks,x,y);
G2=TractionGreenTensor_2D(omega,kp,ks,xx,y);


G(1,:)=G1(1,:)-G2(1,:);

    

