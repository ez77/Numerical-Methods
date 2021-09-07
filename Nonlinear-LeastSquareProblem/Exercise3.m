%%
clear;
% Mesh definition
z = 0:0.001:2;
% Data generation
h = cos(pi.*z)-2*cos(5*pi.*z)+cos(6*pi.*z);
%gaussian noise added
b = h+0.5*randn(size(z));
%fitting curve function
f = @(x) x(1)*cos(x(2)*pi.*z)+x(3)*cos(x(4)*pi.*z)...
    +x(5)*cos(x(6)*pi.*z)-b;
%% Solution using lsqnonlin
%Random initial guess
x0 = [1,1,1,1,-1,5]';%x0=[1,1,1,1,1,5]';x0=[1,1,1,1,-1,5]';
xsol=lsqnonlin(f,x0);
%Plot of the data fitting the curve
plot(z,h,'ko',z,xsol(1)*cos(xsol(2)*pi.*z)+...
    xsol(3)*cos(xsol(4)*pi.*z)...
    +xsol(5)*cos(xsol(6)*pi.*z),'b-');
hold on;
%% Solving using levenberg-marquardt alg.
options.Algorithm = 'levenberg-marquardt';
xsol2 = lsqnonlin(f,x0,[],[],options);
%Plot
plot(z,xsol2(1)*cos(xsol2(2)*pi.*z)+...
    xsol2(3)*cos(xsol2(4)*pi.*z)...
    +xsol2(5)*cos(xsol2(6)*pi.*z),'r.');
%% Solving using fsolve
%Multiplication between the jacobian transpose*the function
%f(x)
g = @(x)[cos(x(2)*pi.*z);-x(1)*pi.*z.*sin(x(2)*pi.*z);...
    cos(x(4)*pi.*z);-x(3)*pi.*z.*sin(x(4)*pi.*z);...
    cos(x(6)*pi.*z);-x(5)*pi.*z.*sin(x(6)*pi.*z)]...
    *(x(1)*cos(x(2)*pi.*z)+x(3)*cos(x(4)*pi.*z)...
    +x(5)*cos(x(6)*pi.*z)-b)';
x02 = [1,1,1,1,-1,5]';
options = optimset('TolX',10^-8);
xsol3 = fsolve(g,x02,options);
%check if the g(x) evaluated ==0
solf = [cos(xsol3(2)*pi.*z);-xsol3(1)*pi.*z.*sin(xsol3(2)*pi.*z);...
    cos(xsol3(4)*pi.*z);-xsol3(3)*pi.*z.*sin(xsol3(4)*pi.*z);...
    cos(xsol3(6)*pi.*z);-xsol3(5)*pi.*z.*sin(xsol3(6)*pi.*z)]...
    *(xsol3(1)*cos(xsol3(2)*pi.*z)+xsol3(3)*cos(xsol3(4)*pi.*z)...
    +xsol3(5)*cos(xsol3(6)*pi.*z)-b)';
%Plot
plot(z,xsol3(1)*cos(xsol3(2)*pi.*z)+...
    xsol3(3)*cos(xsol3(4)*pi.*z)...
    +xsol3(5)*cos(xsol3(6)*pi.*z),'m.-');
%% Solving by Newton method with fin.diff. const. jacob.
n = length(z);
E = eye(6,n);
x03 = [1,1,1,1,-1,5]';
%Newton method with constant jacobian
gNew = @(x)[cos(x(2)*pi.*z);-x(1)*pi.*z.*sin(x(2)*pi.*z);...
    cos(x(4)*pi.*z);-x(3)*pi.*z.*sin(x(4)*pi.*z);...
    cos(x(6)*pi.*z);-x(5)*pi.*z.*sin(x(6)*pi.*z)]*(x(1)*cos(x(2)*pi.*z)+x(3)*cos(x(4)*pi.*z)...
    +x(5)*cos(x(6)*pi.*z)-b)';
m = 6;
Jg = zeros(m);
Eg = eye(m);
delta = 0.001;
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
plot(z,xsolNew(1)*cos(xsolNew(2)*pi.*z)+...
    xsolNew(3)*cos(xsolNew(4)*pi.*z)...
    +xsolNew(5)*cos(xsolNew(6)*pi.*z),'g.-');
hold off;