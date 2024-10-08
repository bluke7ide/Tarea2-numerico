
%% Inciso (b)

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

%% Inciso (c)

% Valores
c = [-1/2, 3/4, -3/2, 1/4];
x0 = 1;

% Verificación con el algoritmo 
disp(evalCheb(c, x0)); 

%% Inciso (d)

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

