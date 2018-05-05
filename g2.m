function y=g2(omega,n)
for i=1:n-1
    y(i)=2*cos((n-i)*(omega-0.02*(omega-pi)^3));
end
y(n)=1;
end