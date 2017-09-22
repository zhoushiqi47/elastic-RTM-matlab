% 在频域下检验parseval等式
n_recv = 1201;
lamda=1/2;

mu=1/4;
omega=pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);
w2 = 1/omega^2 ;
n=size(x);
G=zeros(4,n(2));

l=10*(kp+ks);

kp=kp*(1+0.0000001i);
ks=ks*(1+0.0000001i);
receiver = zeros(2,n_recv);
receiver(1,:) = linspace(-a,a,n_recv);
kk=50;

% 设置间接变量

mus=@(t) (ks^2-t.^2).^(1/2);
mup=@(t) (kp^2-t.^2).^(1/2);
beta=@(t) ks^2-2*t.^2;
delta=@(t) beta(t).^2 + 4*t.^2.*mus(t).*mup(t);
forec1=@(t) 1i./(2*pi.*mu.*delta(t));
forec2=w2*mu/2;

Esy=@(t) exp(1i*mus(t)*10);
Epy=@(t) exp(1i*mup(t)*10);
Esx=@(t) exp(1i*t*receiver(1,kk));

As11=@(t) forec1(t).*mus(t).*beta(t).*Esy(t);
As21=@(t) forec1(t).*2.*mus(t).*mup(t).*t.*Esy(t);

Bs11=@(t) forec2*beta(t).*Esy(t);
Bs21=@(t) forec2*2*mus(t).*t.*Esy(t);


C11=@(t) Bs11(t).*conj(As11(t));



G1=integral(@(t) C11(t) ,-l,l);
Gs1=integral(@(t) As11(t).*Esx(t) ,-l,l);
Trs1=integral(@(t) Bs11(t).*Esx(t)./(2*pi) ,-l,l);
ans1=Gs1-Gs(1,kk)
ans2=Trs1-grs(1,kk)


