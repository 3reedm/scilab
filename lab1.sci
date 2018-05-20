// Решение уравнения теплопроводности методом разностных схем
//       ∂u / ∂t - a ^ 2 * ∂ ^ 2 u / ∂t ^ 2 = f(x, t) 
// a ^ 2 = λ / cρ, 0 < x < L, 0 < t < T, (T - время, любое число до ∞)
// u(0, t) = myu(t) = const, u(L, t) = nyu(t) = const, u(x, 0) = fi(x), 0 < x < L,
// a ^ 2 - коэффициент температуропроводности,
// λ - коэффициент теплопроводности материала стержня,
// c - удельная теплоемкость, ρ - плотность.

// Правая часть дифференциального уравнения.
function y = f(x, t)
    y = 0;
endfunction

// Начальное условие
function y = fi(x)
    y = 2 * sin(2 * %pi * x);
endfunction

// Условие на левой границе
function y = myu(t)
    y = 0;
endfunction

// Условие на правой границе
function y = nyu(t)
    y = 0;
endfunction

// Уравнение метода Фурье u(x, t) = sum(n = 1, ∞) {2 / L * integral(0, L){fi(s) * sin(pi * n * s / L)ds} * sin(pi * n * x / L) * exp (-(a^2) * ((pi * n / L) ^ 2) * t)}
function [uaij] = four(N, K, L, a, x, t)
    h = L / N;
    for i = 1 : N + 1
        X(i) = (i - 1) * h;
    end
    
    uaij = 0;
    for i = 1 : 100
        integ = 0;
        
        for j = 1 : N + 1 
            integ = integ + fi(X(j)) * sin(%pi * i * X(j) / L);
        end
        
        uaijn = 2 / L * integ * sin(%pi * i * x / L) * exp(-(t * (%pi * i * a / L) ^ 2));
   
        uaij = uaij + uaijn;  
    end
endfunction

// Функция решения параболического уравнения методом 
// сеток с помощью явной разностной схемы.
// N - количество участков, на которые разбивается интервал по х (0, L); 
// K - количество участков, на которые разбивается интервал по t (0, T); 
// a - параметр дифференциального уравнения теплопроводности.
// n - число элементов ряда для точного решения через метод Фурье
// Функция возвращает матрицу решений u, матрицу решений ua (точное решение) и вектора x, t
function [u, ua, x, t] = parabol(N, K, L, T, a)
    // Вычисляем шаг по х
    h = L / N;

    // Вычисляем шаг по t
    delta = T / K;

    // Формируем массив x и первый столбец матрицы 
    // решений U из начального условия
    for i = 1 : N + 1
        x(i) = (i - 1) * h;
        u(i, 1) = fi(x(i));
    end
    
    // Формируем массив t, первую и последнюю строку
    // матрицы решений U из граничных условий
    for j = 1 : K + 1
        t(j) = (j - 1) * delta;
        u(1, j) = myu(t(j));
        u(N + 1, j) = nyu(t(j));
    end
    
    // Формируем матрицу решений u с помощью явной разностной схемы
    gam = (a ^ 2) * delta / (h ^ 2);
    for j = 1 : K
        for i = 2 : N
            u(i, j + 1) = gam * u(i - 1, j) + (1 - 2 * gam) * u(i, j) + gam * u(i + 1, j) + delta * f(x(i), t(j));
        end
    end
    
    // Формируем матрицу решений ua с помощью метода фурье
    for j = 1 : K + 1
        for i = 1 : N + 1
            ua(i, j) = four(N, K, L, a, x(i), t(j));
        end
    end
endfunction

// Вызов функции
[U, UA, X, T] = parabol(50, 50, 2, 2, 1);

// Нормирование U
U = U / max(U) * max(UA); 

// Построение графиков
mtlb_hold('on');
surf(X, T, U');
surf(X, T, UA');
title('PARABOLIC EQUATION (TWO VARIANTS)');
xlabel('X');
ylabel('T');
funcprot(0);
