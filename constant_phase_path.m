figure;ezplot('sqrt((sqrt((1-x*x+y*y)^2+4*x*x*y*y)+(1-x*x+y*y))/2)*sin(5*pi/4)+x*cos(5*pi/4)-1',[-10,10])
hold on
ezplot('-sqrt((sqrt((1-x*x+y*y)^2+4*x*x*y*y)+(1-x*x+y*y))/2)*sin(5*pi/4)+x*cos(5*pi/4)-1',[-10,10])
grid
grid on