%% Forward Algorithm function
% This function returns the vector of solutions for a system 
% Lx = b where L is suposed to be a lower triangular matrix
% with no zeros on the diagonal and is taken as an initial 
% parameter
function sol = ForwardEM(m,b)
    L = m;
    sol = zeros(length(b),1);
    sum = 0; 
    if eig(L)~=0
        sol(1) = b(1)/L(1,1);
        for i = 2:length(b)
            for j = 1:i-1
                sum = L(i,j)*sol(j)+sum;
            end
            sol(i)=(b(i)-sum)/L(i,i);
            sum = 0;
        end
    else
        disp('error')
    end
end