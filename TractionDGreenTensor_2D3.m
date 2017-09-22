function G=TractionDGreenTensor_2D3(omega,kp,ks,y,x)
%% n=(0,1)
%% operator T to x
%% compute the traction green tensor of elastic wave operator in 2D
%% x2=0
%% pressure wave Green Tensor




% 向量形式积分
%y is point source
% 2016 12 16 
tic
n=size(x);
G=zeros(4,n(2));

l=10*(kp+ks);
k1=kp*(1+0.0000001i);
k2=ks*(1+0.0000001i);

% 设置间接变量

Gc1=quadv(@(t) A1(t,x(:,:),y(:,:),k1,k2) ,-l,l );
Gc2=integral(@(t) A2(t,x(:,:),y(:,:),k1,k2) ,-l,l, 'ArrayValued',true);
Gc3=integral(@(t) A3(t,x(:,:),y(:,:),k1,k2) ,-l,l, 'ArrayValued',true);
Gc4=quadv(@(t) A4(t,x(:,:),y(:,:),k1,k2) ,-l,l );

G(1,:)=Gc1;
G(2,:)=Gc2;
G(3,:)=Gc3;
G(4,:)=Gc4;
toc
end

function f=A1(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=mus*mup;
C=t*t;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(B*Es+C*Ep)/(2*pi*delta);
end

function f=A2(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=t*mus;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(B*Es-B*Ep)/(2*pi*delta);
end

function f=A3(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=t*mup;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(B*Es-B*Ep)/(2*pi*delta);
end

function f=A4(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=mus*mup;
C=t*t;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(C*Es+B*Ep)/(2*pi*delta);
end


