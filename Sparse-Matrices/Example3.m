clear all;
n = [100 500 2000];
i=1;
sd = n(i)/2;
A = spdiags([-1*diag(eye(n(i))) 100*diag(eye(n(i))) ...
    -1*diag(eye(n(i)))],[-sd 0 sd],n(i),n(i));
A(2:n(i),1)=1;%All rows 1st column
A(1,2:n(i))=1;%Allcolumns 1st row
A(n(i),1:n(i)-1)=1;%All columns last row
A(1:n(i)-1,n(i))=1;%All rows last column
if isequal(A,A')&& min(eigs((A+A')/2)>0)
    smtric=1;%the matrix is symmetric and positive def.
elseif  min(eigs((A+A')/2)<=0)
    smtric=2;%the matrix is not symmetric nor positive def.
else
    smtric=0;%the matrix is not symmetric but positive def.
end
b = ones(n(i),1);
b(2:2:n(i))=-1;
%% a) Solving by chol in full and sparse mode
A1full = full(A);
R1full = chol(A1full);
R1sp = chol(A);
ta1 = tic;
ya1 = R1full'\b;
xa1 = R1full\ya1;
eta1 = toc(ta1);
ta2 = tic;
ya2 = R1sp'\b;
xa2 = R1sp\ya2;
eta2 = toc(ta2);
%% d) commands symamd symrcm
