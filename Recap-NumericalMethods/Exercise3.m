%%
clear;
%% Matrix definition
n = 100;
v = 1:n; v = (v.^3)';
A = toeplitz(v);
xex = 2*ones(n,1); xex(2:2:end)=-1;
b = A*xex;
%% a) Solution by LU factorization with pivoting
%solution of the linear system
[L,U,P] = lu(A);
ysola = L\(P*b);
xsola = U\ysola;
spy(P); %% Pivoting was made
%% relerror, a priori, a posteriori upper bounds
relerrora = abs(norm(xsola-xex))/abs(norm(xex));
[kappa,ub1,ub2] = apriori(A,b);
%% b) SVD decomposition
[Up,S,V] = svd(A);
rA = rank(A); 
%build Ahat with reduced rank rank(Ahat)<rank(A)
cont = 0;
err = 1; %random just for enter the while
%calculate de vector of the division between sigular values
%of S
Sh = S;
Sdiacond = diag(S(1,1)./S(:,:));
cond1 = Sdiacond>2*10^6;
%Evaluate the condition on the singular values (this is the
%most important condition since the err is dependent on Sh
Sh(cond1,cond1) = 0;
rAh=rank(Sh);
Ah = Up*Sh*V';
%Condition on the error frobenius
err = norm(A-Ah,'fro')/norm(A,'fro');
while err>10^-6
    rAh = rAh+1;
end
%% Minimum norm least square solution after SVD decomp
xsolb = lsqminnorm(Ah,b);
%relative error
relerrorb = abs(norm(Ah*xsolb-A*xex))/abs(norm(A*xex));
%% function for the calculation of a priori errors
function [kappa, ub1,ub2] = apriori(A,b)
    kappa = cond(A);
    % condition for unperturbed right hand side
    cons = norm(1000*eps*ones(length(A)))/norm(A);
    ub1 = kappa*cons/(1-kappa*cons);
    % condition for unperturbed matrix A
    cons2 = norm(1000*eps*ones(length(b)))/norm(b);
    ub2 = kappa*cons2;
end