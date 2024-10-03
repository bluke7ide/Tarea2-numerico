



% ejercicio 3
function [T, val, vec] = tortog(n)
    T = zeros(n+1);
    % se sabe que b-a = 2 por lo que b0 = 0, entonces 
    phis = cell(1, n);
    phis{1} = @(x) 1;
    phis{2} = @(x) x;
    for i = 3:n+1
        phis{i} = @(x) (phis{i-1}(x).*(2*i+1).*x-i.*phis{i-2}(x))./(i+1);
    end

    norm = zeros(n+1, 1);
    norm(1) = 1;
    for i = 1:n
        norm(i+1) = integral(@(x) phis{i+1}(x).^2, -1, 1);
    end

    T(1,1) = integral(@(x) x.*phis{1}(x).^2, -1, 1);

    for i = 1:n
        beta = (integral(@(x) x.*phis{i+1}(x).^2, -1, 1))/norm(i+1);
        alpha = sqrt(norm(i+1)/norm(i));
        T(i+1,i+1) = beta;
        T(i,i+1) = alpha;
        T(i+1,i) = T(i,i+1);
    end

    
    
    % w = zeros(n+1, 1);
    % for i = 1:n+1
    %     w{i} = 
    % end

    
    [vec, val] = eig(T);
end

[T, val, vec] = tortog(2);