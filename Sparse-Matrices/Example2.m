n = [100 500 2000];
i=1;
sd = n(i)/2;
A = spdiags([10*diag(eye(n(i))) -1*diag(eye(n(i)))],[0 sd],n(i),n(i));
A(2:n(i),1)=1;%All rows 1st column
A(1,2:n(i))=1;%Allcolumns 1st row
A(n(i),1:n(i)-1)=1;%All columns last row
A(1:n(i)-1,n(i))=1;%All rows last column
if isequal(A,A')&& min(eigs((A+A')/2)>0)
    smtric=1;
elseif  min(eigs((A+A')/2)<0)
    smtric=2;
else
    smtric=0;
end
b = ones(n(i),1);
b(2:2:n(i))=-1;
%%
t1 = tic;
[Ls,Us]=lu(A);
ys = Ls\b;
xsolsp = Us\ys;
et1 = toc(t1);
%%
t2 = tic;
p = colamd(A);
[Lc,Uc]=lu(A(:,p));
yc = Lc\b(p);
xsolc = Uc\yc;
et2 = toc(t2);
