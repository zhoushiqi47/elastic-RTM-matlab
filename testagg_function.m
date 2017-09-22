x=linspace(0,0.5,100);
y=linspace(0,1,100);
z=zeros(100);
for k=1:100
    for j=1:100
        a=x(k);
        b=y(j);
        z(k,j)=16*(a^4-a^3)*b^3 +(24*a^2-16*a^3)*b^2 -8*a*b +1;
    end
end
figure,imagesc(x,y,z.');colorbar