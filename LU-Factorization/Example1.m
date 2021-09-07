A = toeplitz([10 2 3 4 5]);
xex = [1;1;1;1;1];
b = A*xex;
[L,U]=lu(A);
%% Forward and Backward substitution
soly = ForwardEM(L,b);
solx = BackwardEM(U,soly);
abserror = norm(solx-xex);
%% Ex2
d1 = det(A);
dU = det(U);
errordet = abs(d1-dU);
