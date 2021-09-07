function [Lb,Ub]=bandedLU(A)
    Lb = zeros(size(A,1),size(A,2))...
        +diag(diag(eye(length(A))));
    Ub = zeros(size(A,1),size(A,2))+diag(diag(A,1),1);
    d = diag(A);
    ud = diag(A,1);
    ld = diag(A,-1);
    ub(1) = d(1);
    for i = 2:length(A)
        lb(i-1) = ld(i-1)/ub(i-1);
        ub(i) = d(i)-ud(i-1)*lb(i-1);
    end
    Ub = Ub+diag(diag(eye(size(A)).*ub'));
    Lb = Lb+diag(diag(eye(size(A)-1).*lb'),-1);
end