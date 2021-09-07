%%
clear;
%% Matrix Definition
% Size of the matrix
N = 1300;                     e = ones(N,1);
% Matrix construction
A = spdiags([e e e 10*e e e e],[-25 -7 -2 0 2 7 25],N,N);
% Linear system definition
xex = -1*ones(N,1);         xex(1:2:end) = 0;
b = A*sparse(xex);
%extra check on the form of the matrix
figure(1)
spy(A);
%% a) Solution using the command \
ta = tic;
sola = A\b; 
timea = toc(ta);
%% b) LU factorization, checking pivoting
%solution of the linear system
[L,U,P] = lu(A);
%check if pivoting was carried out (easier checking in
%the command window)
pivotflag = isequal(P,eye(N,N));
%solution options if pivoting was made
if pivotflag == 1
    tb1 = tic;
    ysolb = L\b;
    xsolb = U\ysolb;
    timeb1 = toc(tb1);
    timeb2 = 0;
else
    tb2 = tic;
    ysolbpiv = L\(P*b);
    xsolbpiv = U\ysolb;
    timeb2 = toc(tb2);
    timeb1 = 0;
end
%% c) Solution by Cholesky method
% check if the system is positive definite
if(isequal(A,A') && min(eig((A+A')/2))>0)
    Rc = chol(A); %L=R and U=R'
    tc = tic;
    ysolc = (Rc')\b;
    xsolc = Rc\ysolc; 
    timec = toc(tc);
    % the solution has an increase in the error due to the
    % fact that is a banded matrix therefore some reordering
    % could benefit the method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    spy(Rc,'r-'); %sparsity of the matrix without preordering
    hold on
else
    disp('The matrix is not positive definite');
end
%% d) Solution by Cholesky method with permutation using 
% the command symrcm
% Obtain the permutation indexes that reorders the matrix A
r = symrcm(A);
% Solve by cholesky
% If the permutation is made also the b vector must be
% permutated, at the end the solution will be in the order
% xsold=xex(r)
Rd = chol(A(r,r));
bd = b(r);
spy(Rd,'b-'); %sparcity of matrix after ordering -> is reduced
hold off
td = tic;
ysold = (Rd')\bd;
xsold = Rd\ysold;
timed = toc(td);
%% Total results of time, errors
Totime = [timea;timeb1;timeb2;timec;timed];
%% A posteriori error

