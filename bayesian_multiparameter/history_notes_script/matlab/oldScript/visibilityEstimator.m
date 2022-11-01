theta = 1;
V = 0.9;

nu = 100;

p0 = (1+V*cos(theta))/2;
p1 = (1+V*sin(theta))/2;

a0 = binornd(nu, p0);
a1 = binornd(nu, p1);

f0 = a0/nu;
f1 = a1/nu;

vis = (nu*((2*f0-1)^2+(2*f1-1)^2)-1)/(nu-1)
disp(V^2);