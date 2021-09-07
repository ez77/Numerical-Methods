%% Data fitting 
x = [0 0 .1 .1 0 .1 .2 .2 .2]';
y = [0 .1 0 .1 .2 .2 .2 .1 0]';
x2 = (0:.1:.8)';
y2 = ones(9,1)*.2;
A = [x.^2./10 y.^2./10 -1*ones(9,1)];
A2 = [x2.^2./10 y2.^2./10 -1*ones(9,1)];
b = sum(A,2);
b2 = sum(A2,2);
rA = rank(A);
%% Solution with \
xsol1 = A\b;
xsol1E2 = A2\b2;
%% Solution with normal equations chol method
An = A'*A;
bn = A'*b;
Rc = chol(An);
ysol2 = Rc\bn;
xsol2 = Rc'\ysol2;
%% Solution with QR factorization
[Q,R]=qr(A);
xsol3 = R\(Q'*b);
%% Solution with SVD
[U,S,V]=svd(A);
xsol4 = (U*S*V')\b;
%% Solution with pinv
pInv = pinv(A);
xsol5 = pInv*b;

