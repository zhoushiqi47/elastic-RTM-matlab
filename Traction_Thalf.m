function G=Traction_Thalf(omega,kp,ks,x,y)
%% traction tensor of displacement free green tensor
tic
%x is point source
w2 = 1/ks^2 ;
k1=kp*(1+0.00000001i);
k2=ks*(1+0.00000001i);
n=size(x);
G=zeros(4,n(2));
l=10*(kp+ks);
xx=x;
xx(2,:)=-x(2,:);
G1=TractionGreenTensor_2D(omega,kp,ks,x,y);
G2=TractionGreenTensor_2D(omega,kp,ks,xx,y);
for j=1:n(2)

    G(1,j)=integral(@(t) firstTT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(1,j)+G2(1,j);
    G(2,j)=integral(@(t) secondTT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(2,j)+G2(2,j);
    G(3,j)=integral(@(t) thirdTT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(3,j)-G2(3,j);
    G(4,j)=integral(@(t) forthTT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(4,j)-G2(4,j);

    
end
toc
end
%%

function f=firstTT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=-2i*beta.*exp(Ss*y(2))+2i*beta.*exp(Sp*y(2));
B=2i*Ss.*Sp.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))+1i*beta.*t.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B;
end

function f=secondTT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=4*Ss.*t.*t.*exp(Ss*y(2))-beta.*beta./Sp.*exp(Sp*y(2));
B=2i*Ss.*Sp.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))+1i*beta.*t.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B;
end

function f=thirdTT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=-2i*beta.*exp(Ss*y(2))+2i*beta.*exp(Sp*y(2));
B=-2*t.^2.*Sp.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))-beta.*Sp.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B;
end

function f=forthTT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=4*Ss.*t.*t.*exp(Ss*y(2))-beta.*beta./Sp.*exp(Sp*y(2));
B=-2*t.^2.*Sp.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))-beta.*Sp.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B;
end