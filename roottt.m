function [x,y,y1,y2,xx,tt]=roottt(t)
p=[16*t-16, 24-16*t, -8, 1];

x=roots(p);
y=(1-2.*x).^2 +4*x.*(1-x).^(1/2).*(t-x).^(1/2);
y1=(1-2.*x).^2  ;
y2=4*x.*(1-x).^(1/2).*(t-x).^(1/2) ;
xx=x.^(1/2);
tt=sqrt(t);

