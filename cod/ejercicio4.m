
%% Inciso (f)

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

%% Inciso (g)
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


