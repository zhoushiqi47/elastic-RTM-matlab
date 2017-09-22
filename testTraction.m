function f=testTraction(t,x,y,w2,kp,ks)

Ss=(t.^2-ks.^2).^(1/2);
Sp=(t.^2-kp.^2).^(1/2);
beta=ks.^2-2.*t.^2;

B=2i*Ss.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))+1i*beta./Sp.*t.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
B=2i*Ss.*Sp.*t.*exp(Ss*x(2)+1i*t*(y(1)-x(1)))+1i*beta.*t.*exp(Sp*x(2)+1i*t*(y(1)-x(1)));
f=w2./(4*pi)*B./Sp;
end
