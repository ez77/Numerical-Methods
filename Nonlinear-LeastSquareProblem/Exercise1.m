%% Function definition
clear;
Lx = 2;
Ly = 2;
f = @(x,y) 15+x+2*y-(80/(Lx.^2)).*((x-Lx/2).^2+(y-Ly/2).^2)...
    +300/Ly.^4.*((x-Lx/2).^4+(y-Ly/2).^4);
x0 = [0.5;0.5];
tol = 10^-12;
gamma = 0.001;
%% Plotting function
[X,Y] = meshgrid(0:0.01:2);
figure (1)
surf(X,Y,f(X,Y));
figure (2)
contour(X,Y,f(X,Y));
% From the figures we can say that the minimum is approx at
% [0.22,0.22] and f = 5.2
%% Steepest descent
%Convergence check with hessian matrix
syms x y;
fsym = 15+x+2*y-(80/(Lx^2))*((x-Lx/2)^2+(y-Ly/2)^2)...
    +300/Ly^4*((x-Lx/2)^4+(y-Ly/2)^4);
hmat = hessian(fsym,[x,y]);
values = eig(subs(hmat,[x y],[x0(1) x0(2)])); % the values ar positive so it must converge
values = double(values);
%Calculation of the minimum value
fpx = @(x) 1-(80/Lx.^2).*(2.*(x-Lx/2))...
    +(300/Ly.^4).*(4.*(x-Lx/2).^3);
fpy = @(y) 2-(80/Lx.^2).*(2.*(y-Ly/2))...
    +(300/Ly.^4).*(4.*(y-Ly/2).^3);
gdient = [fpx(x0) fpy(x0)]';
relerror = 1;
iter = 0;
while relerror>tol
    xsol = x0 -gamma*fpx(x0);
    relerror = abs(f(xsol(1),xsol(2))-f(x0(1),x0(2)))...
        /abs(f(x0(1),x0(2)));
    x0 = xsol;
    iter = iter+1;
end
%% Modified gradient method
% Initial parameters definition
x02 = [0.5;0.5];
relerror2 = 1;
iter2 = 0; %#iterations for the mod. grad. method
iterg = 0; %#iterations for the condition on m 
m = 8;  %the value of m must be choosen in order to comply with 
        %Machine error 
while relerror2>=tol 
    gmod = gamma/2^m;
    xsol2 = x02 -gmod*fpx(x02);
    %Condition for m
    upplim = -(gamma/2^(m+1))*norm(fpx(x02))^2;
    cond = f(xsol2(1),xsol2(2))-f(x02(1),x02(2));
    while iterg<1000 && upplim<cond
        if (1/2^m<100*eps)
            m = m-1;
        end
        iterg = iterg +1;
    end
    %Calculation of the relative error
    relerror2 = abs(f(xsol2(1),xsol2(2))-f(x02(1),x02(2)))...
        /abs(f(x02(1),x02(2)));
    x02 = xsol2;
    iter2 = iter2+1;
end
%% Newton method
%Initial parameters definition
x03 = [0.2;0.2];
relerror3 = 1;
iter3 = 0;
while relerror3>=tol
    %deltax = invhessianmatrix*gradient of the function to min.
    deltax = -1*double(subs(hmat,[x y],[x03(1) x03(2)]))\fpx(x03);
    xsol3 = x03 + deltax;
    relerror3 = abs(f(xsol3(1),xsol3(2))-f(x03(1),x03(2)))...
        /abs(f(x03(1),x03(2)));
    x03 = xsol3;
    iter3 = iter3+1;
end
