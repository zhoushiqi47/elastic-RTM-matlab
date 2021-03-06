function G=Elastic_GreenTensor_Thalf2(omega,kp,ks,x,y)
%没有表面波部分
%x is point source
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
for j=1:n(2)

    G(1,j)=G1(1,j)+G2(1,j);
    G(2,j)=G1(2,j)+G2(2,j);
    G(3,j)=G1(3,j)-G2(3,j);
    G(4,j)=G1(4,j)-G2(4,j);

    
end
end

function f=firstG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=2i*Ss.*t.*exp(Ss*y(2))+1i*beta./Sp.*t.*exp(Sp*y(2));
B=2i*Ss.*Sp.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))+1i*beta.*t.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B;
end

function f=secondG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=2*t.^2.*exp(Ss*y(2))+beta.*exp(Sp*y(2));
B=2i*Ss.*Sp.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))+1i*beta.*t.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B;
end

function f=thirdG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=2i*Ss.*Sp.*t.*exp(Ss*y(2))+1i*beta.*t.*exp(Sp*y(2));
B=-2*t.^2.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))-beta.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B;
end

function f=forthG(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;
delta=beta.^2-4*t.^2.*Ss.*Sp;
A=2*t.^2.*exp(Ss*y(2))+beta.*exp(Sp*y(2));
B=-2*t.^2.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))-beta.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(2*pi*delta).*A.*B.*Sp;
end


