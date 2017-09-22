function G=Traction_Dhalf(omega,kp,ks,x,y)
%% traction tensor of displacement free green tensor

%x is point source
w2 = 1/(ks^2) ;
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

    G(1,j)=integral(@(t) firstT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(1,j)-G2(1,j);
    G(2,j)=integral(@(t) secondT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(2,j)-G2(2,j);
    G(3,j)=integral(@(t) thirdT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+integral(@(t) thirdT2(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(3,j)+G2(3,j);
    G(4,j)=integral(@(t) forthT(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(4,j)+G2(4,j);

    
end


end
%%

function f=firstT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
alpha=t.^2-Ss.*Sp;
A=1i*beta.*t./Ss.*exp(Ss*y(2)) + 2i*Sp.*t.*exp(Sp*y(2));
B=1i*t.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))-1i*t.*Ss.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=secondT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
alpha=t.^2-Ss.*Sp;
A=-2*t.*t.*exp(Ss*y(2)) -beta.*exp(Sp*y(2));
B=1i*t.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))-1i*t.*Ss.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=thirdT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
alpha=t.^2-Ss.*Sp;
A=1i*beta.*t./Ss.*exp(Ss*y(2))+2i*Sp.*t.*exp(Sp*y(2)) ;
B=-Sp.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))+t.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=thirdT1(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
alpha=t.^2-Ss.*Sp;
A=1i*beta.*t./Ss.*exp(Ss*y(2)) ;
B=-Sp.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))+t.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=thirdT2(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
alpha=t.^2-Ss.*Sp;
A= 2i*Sp.*t.*exp(Sp*y(2));
B=-Sp.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))+t.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=forthT(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
alpha=t.^2-Ss.*Sp;
A=-2*t.*t.*exp(Ss*y(2)) -beta.*exp(Sp*y(2));
B=-Sp.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))+t.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end
