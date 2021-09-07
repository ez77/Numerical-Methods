clear all;
n = 5000;
A = spdiags([-1*diag(eye(n)) -1*diag(eye(n)) -1*diag(eye(n))...
     7*diag(eye(n)) 1*diag(eye(n)) 1*diag(eye(n))... 
     1*diag(eye(n))],[-3000 -80 -30 0 20 50 2000],n,n);
Af = full(A);
b = ones(n,1);
b(2:2:n)=0;
%% a) solve by \, lu and symrcm for sparse matrix
x1 = Af\b;
[L,U]=lu(Af);
y = L\b;
x2 = U\y;
% Put attention in the order of the solution since for this
% case it was necessary to reorder b but for colamd no. 
r = symrcm(A);
[Ls,Us]=lu(A(r,r));
ys = Ls\b(r);
xs = Us\ys;
%% Error calculation
abserrinfFULL = abs(norm(b,inf)-norm(A*x1,inf));
abserrinfSPARSE = abs(norm(b(r),inf)-norm(A(r,r)*xs,inf));
abserrl2SFULL = abs(norm(b)-norm(A*x1));
abserrl2SPARSE = abs(norm(b(r))-norm(A(r,r)*xs));
