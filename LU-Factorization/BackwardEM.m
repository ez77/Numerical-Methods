%% Backward substitution algorithm
% This functions returns a vector of solutions for the problem Ux = b where
% U is suposed to be an upper triangular matrix taken as an inital
% parameter
function sol = BackwardEM(U,b)
    sol = zeros(length(b),1);
    sum = 0; 
    if eig(U)~=0
        sol(end) = b(end)/U(end,end);
        for i = length(b)-1:-1:1
            for j = i+1:length(b)
                sum = U(i,j)*sol(j)+sum;
            end
            sol(i)=(b(i)-sum)/U(i,i);
            sum = 0;
        end
    else
        disp('error')
    end
end