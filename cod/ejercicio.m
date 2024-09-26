



% ejercicio 3
function [T, val, vec, w, nod] = tortog(n)
    T = zeros(n+1);
    % se sabe que b-a = 2 por lo que b0 = 0, entonces 
    % phis = zeros(n+1,1);
    % phis(1) = 1;
    % phis(2) = x;
    % for i = 3:n+1
    % 
    %     phis(i) = 
    % end

    % for i = 1:n
    %     for j = 1:n
    %        T(i,j) = 
    %        T(i,j+1) = 
    %        T(i+1,j) = T(i,j+1);
    %     end
    % end
    % 
    % T(n+1, n+1) = 



    
    [vec, val] = eig(T);
end

tortog(0)