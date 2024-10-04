


%%
% ejercicio 3

% parte i
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
        T(k,k+1) = sqrt(k/(2*k + 1));
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

% parte j

% [~, nod, w] = tortog(0);
% [~, nod, w] = tortog(1);

% parte k
n = 5;
[T, nod, w] = tortog(n-1);
f = @(x) x^8 + 2*x^2 + x;
int = 0;
for k = 1:n
    int = int + w(k)*f(nod(k));
end

intf = @(x) x^9/9 + 2/3*x^3 + 0.5*x^2;
explicit = intf(1) - intf(-1);

%%