



% ejercicio 3
function [T, nod, w] = tortog(n)
    T = zeros(n+1);
    for k = 1:n
        T(k,k+1) = sqrt(k/(k+1));
        T(k+1,k) = T(k,k+1);
    end
    
    % w = zeros(n+1, 1);
    % for i = 1:n+1
    %     w{i} = 
    % end

    
    [w, nod] = eig(T);
end

[T, nod, w] = tortog(1);