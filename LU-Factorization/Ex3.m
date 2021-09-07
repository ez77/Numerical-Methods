%% Matrix definitions
B = [10 3 4 1 -1;
    1 30 3 4 5;
    1 2 50 4 5;
    1 2 3 30 5;
    1 2 3 4 10];
A = kron(eye(200),B);
xex = diag(eye(length(A)));
bmat = A*xex;
%% lu matlab command and functions
[L,U,P]=lu(A);
tic
solmaty = L\bmat;
solmatx = U\solmaty;
toc
soly = ForwardEM(L,bmat);
solx = BackwardEM(U,soly);
toc
%% Error calculation
abserrormat = norm(solmatx-xex);
relerrormat = norm(solmatx-xex)/norm(xex);
abserror = norm(solx-xex);
relerror = abserror/norm(xex);