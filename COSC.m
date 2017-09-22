function f=COSC(x)
a=real(x);
b=imag(x);
f=1i*(exp(-b)-exp(b))./2.*sin(a);