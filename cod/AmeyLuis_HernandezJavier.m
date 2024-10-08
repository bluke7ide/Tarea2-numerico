%% ejercicio 3

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
disp(int)

intf = @(x) x^9/9 + 2/3*x^3 + 0.5*x^2;
disp(intf(1) - intf(-1))

% parte l

function c = gausscomp(a, b, m, n, f)
    %
    %   Función que calcula la cuadratura de Gauss compuesta, por
    %   medio de subintervalos m, con un nivel n, en el intervalo
    %   [a,b] de la función f
    %
    %   Inputs
    %       a - intervalo izquierdo
    %       b - intervalo derecho
    %       m - la cantidad de subintervalos
    %       n - el nivel de las cuadraturas de Gauss
    %       f - la función a integrar
    % 
    %   Outputs
    %       c - el valor de la integral

    % Realizamos el ancho de banda y los espaciamientos
    h = (b - a) / m;  
    xj = linspace(a, b, m+1);  
    
    % Hacemos la doble sumatoria de la cuadratura compuesta
    [~, nod, w] = tortog(n);  
    c = 0;
    for j = 1:m
        for k = 1:n+1
            c = c + w(k) * f( (xj(j) + xj(j+1))/ 2 + (h / 2)*nod(k) );
        end
    end
    c = c * h / 2;
end

% a modo de comprobación
disp(gausscomp(-1, 1, 1, 5, f))

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

%% ejercicio 4
% Inciso (f)

function D = calcularMatrizD(n)
    %
    %   Función que calcula la matriz D con las derivadas del polinomio de
    %   de Lagrange evaluadas en los nodos de chebyshev.
    %
    %   Inputs:
    %       n - grado de la función.
    %   
    %   Outputs:
    %       D - matriz con las derivadas del polinomio de Lagrange
    %           evaluadas sobre los nodos.
    % 

    % Calcula los nodos de Chebyshev
    x_chebyshev = cos((0:n)' * pi / n); 

    % Calcula los lambdas para todo x nodo de Chebyshev
    lambdas = zeros(n+1, 1);

    for k = 1:n+1
        prod = 1;
        for i = 1:n+1
            if i ~= k
                prod = prod * (x_chebyshev(k) - x_chebyshev(i));
            end
        end
        lambdas(k) = 1 / prod;
    end

    % Inicializa la matriz D
    D = zeros(n+1, n+1); 

    % Llena la matriz D en las entradas i!=k
    for i = 1:n+1
        for k = 1:n+1
            if i ~= k
                % Caso i != k
                D(i, k) = lambdas(k) / (lambdas(i) * (x_chebyshev(i) - x_chebyshev(k)));
            end
        end
    end

    % Llena la matriz D en las entradas i=k
    for i = 1:n+1
        D(i, i) = -sum(D(i, :));  
    end

end

% Ejecución de la función
n = 20;
f = @(x) 1 ./ (1 + 16*x.^2);

% Nodos de Chebyshev
x_chebyshev = cos((0:n)' * pi / n);

% Evaluación de f en los nodos
fx = f(x_chebyshev);

% Calcula la matriz D
D = calcularMatrizD(n);

% Derivada exacta de f
df_exacta = @(x) -32 * x ./ (1 + 16*x.^2).^2;
x_continuo = linspace(-1,1);

% Derivada aproximada
df_aprox = D * fx;

% Gráfica
figure;
plot(x_continuo, df_exacta(x_continuo), 'r', 'DisplayName', 'Derivada exacta', 'LineWidth', 2);
hold on;
plot(x_chebyshev, df_aprox, 'b.', 'DisplayName', 'Derivada aproximada', 'MarkerSize', 15);
legend;
xlabel('x');
ylabel('f''(x)');
title('Comparación entre derivada exacta y aproximada');
grid on;

% Cuantificación del error
error_df = abs(df_exacta(x_chebyshev) - df_aprox);
disp(error_df);


% Inciso (g)
n = 1000;
errores = zeros(1, n);

for k = 1:n
    % Calcula la matriz D para cada k
    D = calcularMatrizD(k);
    x_chebyshev = cos((0:k)' * pi / k);
    
    % Evalúa la función y la derivada
    fx = f(x_chebyshev);
    df_aprox = D * fx;
    dfx_exacta = df_exacta(x_chebyshev);
    
    % Norma infinito
    errores(k) = max(abs(df_aprox - dfx_exacta));
end

% Graficar el error
figure;
semilogy(1:n, errores(1:n), '.-');
xlabel('n');
ylabel('Error');
title('Error en la derivada aproximada en función de n');
grid on;

% Error mínimo
[val_min, pos_min] = min(errores);
disp(pos_min);
disp(val_min);


%% ejercicio 5
% Inciso (b)

function y0 = evalCheb(c, x0)
    %
    %   Función que calcula el valor del polinomio evaluado en x0 de manera
    %   recursiva definida por el Algoritmo 1, mediante los coeficientes
    %   de Chebyshev.
    %
    %   Inputs:
    %       c - vector de coeficientes de Chebyshev
    %       x0 - punto donde se evalúa el polinomio p(x)
    %
    %   Outputs:
    %       y0 - valor del polinomio en x0 según el Algoritmo 1.
    % 

    % Valores iniciales
    n = length(c)-1;  
    b_n2 = 0;  
    b_n1 = 0;  

    % Recurrencia
    for k = n:-1:1
        b_k = c(k+1) + 2 * x0 * b_n1 - b_n2;  
        b_n2 = b_n1;  
        b_n1 = b_k;   
    end

    % Evaluación de p(x0)
    y0 = c(1) + x0 * b_n1 - b_n2;  

end


% Inciso (c)

% Valores
c = [-1/2, 3/4, -3/2, 1/4];
x0 = 1;

% Verificación con el algoritmo 
disp(evalCheb(c, x0)); 

% Inciso (d)

% Valores
n = 10^8;
c = rand(n,1);
x0 = 0.1;

% Número de ejecuciones
num_iteraciones = 10;

% evalCheb
tiempos_cheb = zeros(num_iteraciones, 1);
for i = 1:num_iteraciones
    tic;  
    evalCheb(c, x0);  
    tiempos_cheb(i) = toc;  
end
t_cheb = mean(tiempos_cheb);
disp(t_cheb)

% polyval
tiempos_polyval = zeros(num_iteraciones, 1);
for i = 1:num_iteraciones
    tic;  
    polyval(c, x0);  
    tiempos_polyval(i) = toc;  
end
t_polyval = mean(tiempos_polyval);
disp(t_polyval)