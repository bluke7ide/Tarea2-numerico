
%% Inciso (f)

function D = calcularMatrizD(n)
    % n: Número de nodos

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
        lambdas(k) = prod;
    end

    % Inicializa la matriz D
    D = zeros(n+1, n+1); 

    % Llena la matriz D en las entradas i≠k
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
x = cos((0:n)' * pi / n);

% Cálcula la matriz D
D = calcularMatrizD(n);

% Derivada exacta de f
df_exacta = @(x) -32 * x ./ (1 + 16*x.^2).^2;

% Evaluación de f en los nodos
fx = f(x);

% Derivada aproximada
df_aprox = D * fx;

% Gráfica
figure;
plot(x, df_exacta(x), 'r', 'DisplayName', 'Derivada exacta');
hold on;
plot(x, df_aprox, 'b.', 'DisplayName', 'Derivada aproximada', 'MarkerSize', 15);
legend;
xlabel('x');
ylabel('f''(x)');
title('Comparación entre derivada exacta y aproximada');
grid on;


%% Inciso (g)

errores = zeros(1, n);

for k = 1:n
    % Calcula la matriz D para cada k
    D = calcularMatrizD(k);
    x = cos((0:k)' * pi / k);
    
    % Evalúa la función y la derivada
    fx = f(x);
    df_aprox = D * fx;
    df_exacta_x = df_exacta(x);
    
    % Norma infinito
    errores(k) = max(abs(df_aprox - df_exacta_x));
end

% Graficar el error
figure;
plot(1:n, errores(1:n), '-o');
xlabel('n');
ylabel('Error');
title('Error en la derivada aproximada en función de n');
grid on;


