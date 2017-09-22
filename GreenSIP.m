function G=GreenSIP(k,x,y)
mu=@(t) (k^2-t.^2).^(1/2);
l=10*k;
if x(2,:)>y(2,:)
    G=1/(4*pi)*integral(@(t) 1i./mu(t).*exp(1i*mu(t)*(x(2,:)-y(2,:))+1i*t*(x(1,:)-y(1,:))),-l+1i,l-1i,'Waypoints',[-1+1i,1-1i]);
else
    G=1/(4*pi)*integral(@(t) 1i./mu(t).*exp(1i*mu(t)*(y(2,:)-x(2,:))+1i*t*(x(1,:)-y(1,:))),-l+1i,l-1i,'Waypoints',[-1+1i,1-1i]);
end
