n=1000;
x=linspace(50,100,n);
x=exp(x);
ks=2*pi;
E=@(t) exp((10*i*sqrt(1-t^2)+1i*t*x));
y=ks*integral(@(t) E(t) ,-1,1, 'ArrayValued',true);
z1=log(x);
z2=log(abs(y));
figure;
plot(z1,z2);
t=(z2(n)-z2(1))/(z1(n)-z1(1))
