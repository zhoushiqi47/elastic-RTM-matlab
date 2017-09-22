function G=Elastic_GreenTensor_Dhalf(omega,kp,ks,x,y)
tic
%displacement free Green Tensor
%x is point source
w2 = 1/omega^2 ;
k1=kp*(1+0.00000001i);
k2=ks*(1+0.00000001i);
n=size(x);
G=zeros(4,n(2));
l=10*(kp+ks);
xx=x;
xx(2,:)=-x(2,:);
G1=Elastic_GreenTensor_2D(omega,kp,ks,x,y);
G2=Elastic_GreenTensor_2D(omega,kp,ks,xx,y);
for j=1:n(2)

    G(1,j)=integral(@(t) firstGG(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(1,j)-G2(1,j);
    G(2,j)=integral(@(t) secondGG(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(2,j)-G2(2,j);
    G(3,j)=integral(@(t) thirdGG(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(3,j)+G2(3,j);
    G(4,j)=integral(@(t) forthGG(t,x(:,j),y(:,j),w2,k1,k2) ,-l,l)+G1(4,j)+G2(4,j);

    
end
toc
end

function f=firstGG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);

alpha=t.^2-Ss.*Sp;
A=-1i*t.*exp(Ss*y(2)) + 1i*t.*exp(Sp*y(2));
B=1i*t.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))-1i*t.*Ss.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=secondGG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);

alpha=t.^2-Ss.*Sp;
A=-t.*t./Ss.*exp(Ss*y(2)) + Sp.*exp(Sp*y(2));
B=1i*t.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))-1i*t.*Ss.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=thirdGG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);

alpha=t.^2-Ss.*Sp;
A=-1i*t.*exp(Ss*y(2)) + 1i*t.*exp(Sp*y(2));
B=-Sp.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))+t.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end

function f=forthGG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);

alpha=t.^2-Ss.*Sp;
A=-t.*t./Ss.*exp(Ss*y(2)) + Sp.*exp(Sp*y(2));
B=-Sp.*Ss.*exp(Sp*x(2)+1i*t*(y(1)-x(1)))+t.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*alpha).*A.*B;
end
