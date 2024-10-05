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
        T(k,k+1) = k/sqrt(4*k^2 - 1);
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

[~, nod, w] = tortog(0);
disp(nod)
disp(w)
[~, nod, w] = tortog(1);
disp(nod)
disp(w)

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

% parte l

function c = gausscomp(a, b, m, n, f)
    h = (b - a) / m;  
    xj = linspace(a, b, m+1);  
    
    [~, nod, w] = tortog(n);  
    
    c = 0;
    
    for j = 1:m
        for k = 1:n+1
            c = c + w(k) * f( (xj(j) + xj(j+1))/ 2 + (h / 2)*nod(k) );
        end
    end
    c = c * h / 2;
end

gausscomp(-1, 1, 1, 5, f)

% parte m

g = @(x) sin(100.*pi.*x).*((1-x).^0.5).*log(1-x);
m = zeros(100, 4);
ind = [0, 2, 5, 8];
for j = 1:4
    for i = 1:100
        m(i,j) = abs(gausscomp(0, 1, i, ind(j), g) - 0.000819761237123984);
    end
end

x = 1:100;  
plot(x, m(:, 1), '-', x, m(:, 2), '-', x, m(:, 3), '-', x, m(:, 4), '-');

xlabel('m');
ylabel('error abs');
title('Gráfico de los errores absolutos ');
legend('n = 0', 'n = 2', 'n = 5', 'n = 8 ');
grid on;
