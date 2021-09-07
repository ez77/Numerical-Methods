%% Matrix definition
n = 300;
Ah = 4*diag(diag(eye(n)))-1*diag(diag(eye(n-1)),1)...
    -1*diag(diag(eye(n-1)),-1);
Ah(1,n) = -1/2;             Ah(n,1) = Ah(1,n);
bh = ones(300,1);
A = kron(ones(10,1),Ah);    b = kron((1:10)',bh);
rank(A)
As = sparse(A);             bs = sparse(b);
%% Solution with \
xsol1 = A\b;
%% Solution with normal equations chol method
An = A'*A;          bn = A'*b;
Ans = As'*As;       bns = As'*bs;
Rc = chol(An);      Rcs = chol(Ans);
ysol2 = Rc\bn;      xsol2 = Rc'\ysol2;
ysol2sp = Rcs\bns;  xsol2sp = Rcs'\ysol2sp;
%% Solution with QR factorization
[Q,R]=qr(A);
xsol3 = R\(Q'*b);
%% Solution with SVD
% [U,S,V]=svd(A);
% xsol4 = (U*S*V')\b;
