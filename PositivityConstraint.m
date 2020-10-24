

rho = -0.2;
xi = 3;
hx = 1/1000;
hv = 1/10000;


pi = 0.3:0.01:0.7;
x = 30:1:70;


P = length(pi);
X = length(x);

Y1 = zeros(P,X);
Y2 = zeros(P,X);

for i = 1:1:P
for j = 1:1:X
    Y1(i,j) = pi(j)*x(i)/hx - abs(rho*xi)/hv;
    Y2 (i,j)= xi/hv - abs(x(i)*rho*pi(j))/hx;
end
end

figure()
surf(x,pi,Y1)
hold on
surf(x,pi,Y2)
hold on
surf(x,pi,zeros(P,X))




