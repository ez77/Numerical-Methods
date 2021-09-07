%% 
clear;
tp = 1:0.01:3;
%Data generation
h = tp.^3+tp.^2-1;
%gaussian noise added
b = h+0.5*randn(size(tp));
%fitting curve function
f = @(x) tp.^x(1)+tp.^x(2)+x(3)-b;
%Random initial guess
x0 = [1,1,1]';
xsol=lsqnonlin(f,x0);
%Plot of the data fitting the curve
plot(tp,h,'ko',tp,tp.^xsol(1)+tp.^xsol(2)+xsol(3),'b-');
hold on;
%% Using levenberg-marquardt method
options.Algorithm = 'levenberg-marquardt';
xsol2 = lsqnonlin(f,x0,[],[],options);
plot(tp,tp.^xsol2(1)+tp.^xsol2(2)+xsol2(3),'r.');
%% Using fsolve
g = @(x)[log(tp).*tp.^x(1);log(tp).*tp.^x(2);tp./tp]...
    *(tp.^x(1)+tp.^x(2)+x(3)-b)';
x02 = [2.5,1,1]';
options = optimset('TolX',10^-8);
xsol3 = fsolve(g,x02,options);
%check if the g(x) evaluated ==0
solf = [log(tp).*tp.^xsol3(1);log(tp).*tp.^xsol3(2);tp./tp]...
    *(tp.^xsol3(1)+tp.^xsol3(2)+xsol3(3)-b)';
plot(tp,tp.^xsol3(1)+tp.^xsol3(2)+xsol3(3),'m*');
%% finite diff. constant Jacobian with Newton method
n = length(tp);
Jact = zeros(3,n);
E = eye(3,n);
delta = 0.01;
x03 = [2.8,1,1]';
for i = 1:3
    %transpose jacobian finite difference for constucting
    %g(x)
    Jact(i,:)=(f(x03+E(:,i)*delta)-f(x03))/delta;
end
%Newton method with constant jacobian
gNew = @(x)Jact*(tp.^x(1)+tp.^x(2)+x(3)-b)';
m = 3;
Jg = zeros(m);
Eg = eye(m);
delta = 0.01;
for i = 1:m
    %jacobian finite difference for constucting
    %g(x)
    Jg(:,i)=(gNew(x03+E(:,i)*delta)-gNew(x03))/delta;
end
relerror = 1;
tol = 10^-8;
iter = 0;
while relerror>=tol
    deltax = -1*(Jg\gNew(x03));
    xsolNew = x03+deltax;
    relerror = abs(xsolNew-x03)/abs(x03);
    x03 = xsolNew;
    iter = iter+1;
end
plot(tp,tp.^xsolNew(1)+tp.^xsolNew(2)+xsolNew(3),'g.-');