



% ejercicio 3
function [T, nod, w] = tortog(n)
    %
    %   Función que recibe el tamaño de la matriz T y calcula los nodos
    %   y pesos de la cuadratura de Gauss
    %
    %   Inputs:
    %       n   - el tamaño de la matriz T
    %
    %   Outputs:
    %       T   - matriz de los alpha y beta demostrada
    %       nod - los nodos de la cuadratura de Gauss
    %       w   - los pesos de la cuadratura de Gauss

    % Iniciamos los valores
    T = zeros(n+1);
    nod = zeros(n+1, 1);
    w = zeros(n+1, 1);

    % Rellenamos los alpha
    for k = 1:n
        T(k,k+1) = sqrt(k/(k+1));
        T(k+1,k) = T(k,k+1);
    end

    % Calculamos los valores y vectores propios
    [vec, val] = eig(T);
    
    % Calculamos los nodos y pesos de la cuadratura de Gauss
    for i = 1:n+1
        w(i) = 2*(vec(1, i))^2; 
        nod(i) = val(i,i);
    end
end

[T, nod, w] = tortog(1);