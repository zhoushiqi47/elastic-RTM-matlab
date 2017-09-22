labmda=pi;
mu=pi;
tt=linspace(0.1,2,100);
t=tt-0.1i;
beta=(1-2*t.*t);
mu1=(1-t.*t).^(1/2);
mu2=(mu/(labmda+2*mu)-t.*t).^(1/2);
yt=beta.^2 + 4*mu1.*mu2.*t.*t;
yy=abs(yt);
plot(tt,yy);