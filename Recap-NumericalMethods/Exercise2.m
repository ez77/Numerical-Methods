%%
clear;
%% Problem parameters
n = 5;
xex = 2*diag(eye(n));
x0 = 1.87*diag(eye(n));
%% a) Solving by fsolve
f = @prob;
options = optimset('TolX',10^-12);
xsola = fsolve(f,x0,options);
% The solution diverges for the parameter x0<1.87 approx so
% the ring of convergence is really small
%% b) Newton method with constant Jacobian
%% c) Newton method with forward finite difference
%% Function definitions of g(x) and f(x)
function f = prob(x)
    n = 5;
    A = diag(4*diag(eye(n)))+diag(-1*diag(eye(n-1)),1)...
        +diag(-1*diag(eye(n-1)),-1);
    for i=1:n
        gfun(i) = 2-2.*exp(x(i).^2);
    end
    sum = 0;
    for i = 1:n
        for j = 1:n
            sum =A(i,j)*(x(j)-2)+sum;
        end
        f(i) = gfun(i)*sum;
    end
end