function G=Elastic_GreenTensor_Thalf1(omega,kp,ks,x,y)
% 向量形式积分
%x is point source
% 2016 6 22 

w2 = 1/omega^2 ;
k1=kp*(1+0.000001i);
k2=ks*(1+0.000001i);
n=size(x);
G=zeros(4,n(2));

l=10*(kp+ks);

xx=x;
xx(2,:)=-x(2,:);
G1=Elastic_GreenTensor_2D(omega,kp,ks,x,y);
G2=Elastic_GreenTensor_2D(omega,kp,ks,xx,y);

Ss=@(t) (t.^2-k2.^2).^(1/2);
Sp=@(t) (t.^2-k1.^2).^(1/2);
beta=@(t) k2.^2-2.*t.^2;
delta=@(t) beta(t).^2-4*t.^2.*Ss(t).*Sp(t);

A1=@(t) w2./(2*pi*delta(t)).*(2i*Ss(t).*t.*exp(Ss(t)*y(2,:))+1i*beta(t)./Sp(t).*t.*exp(Sp(t)*y(2,:)));
A2=@(t) w2./(2*pi*delta(t)).*(2*t.^2.*exp(Ss(t)*y(2,:))+beta(t).*exp(Sp(t)*y(2,:)));

B1=@(t) 2i*Ss(t).*Sp(t).*t.*exp(Ss(t)*x(2,:)+1i*t*(y(1,:)-x(1,:)))+1i*beta(t).*t.*exp(Sp(t)*x(2,:)+1i*t*(y(1,:)-x(1,:)));
B2=@(t) -2*t.^2.*Sp(t).*exp(Ss(t)*x(2,:)+1i*t*(y(1,:)-x(1,:)))-beta(t).*Sp(t).*exp(Sp(t)*x(2,:)+1i*t*(y(1,:)-x(1,:)));



    G(1,:)=integral(@(t) A1(t).*B1(t) ,-l,l, 'ArrayValued',true)+G1(1,:)+G2(1,:);
    G(2,:)=integral(@(t) A2(t).*B1(t) ,-l,l, 'ArrayValued',true)+G1(2,:)+G2(2,:);
    G(3,:)=integral(@(t) A1(t).*B2(t) ,-l,l, 'ArrayValued',true)+G1(3,:)-G2(3,:);
    G(4,:)=integral(@(t) A2(t).*B2(t) ,-l,l, 'ArrayValued',true)+G1(4,:)-G2(4,:);
    

end


    

