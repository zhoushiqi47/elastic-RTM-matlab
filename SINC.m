function f=SINC(x)
a=real(x);
b=imag(x);
f=(exp(b)+exp(-b))./2.*sin(a);