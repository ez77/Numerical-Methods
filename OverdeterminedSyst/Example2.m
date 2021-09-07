%% Matrix definition
Ah = 2*diag(diag(eye(20)))-1*diag(diag(eye(19)),1)...
    -1*diag(diag(eye(19)),-1);
bh = ones(20,1);
A = kron(ones(10,1),Ah);
b = kron((1:10)',bh);
As = sparse(A);
bs = sparse(b);
rank(A)
%% Solution with \
xsol1 = A\b;
%% Solution with normal equations chol method
An = A'*A;
bn = A'*b;
Ans = As'*As;
bns = As'*bs;
Rc = chol(An);
Rcs = chol(Ans);
ysol2 = Rc\bn;
xsol2 = Rc'\ysol2;
ysol2sp = Rcs\bns;
xsol2sp = Rcs'\ysol2sp;
%% Solution with QR factorization
[Q,R]=qr(A);
xsol3 = R\(Q'*b);
%% Solution with SVD
[U,S,V]=svd(A);
xsol4 = (U*S*V')\b;
%% Solution with pinv
pInv = pinv(A);
xsol5 = pInv*b;
