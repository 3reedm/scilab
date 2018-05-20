// Сохранение графиков
function[]=save_graphics(b,n)
    a = 0.001:0.001:b-0.001
    l = length(a)
    
    for i=1:l
        [Pi(i), N(i), N_m(i), Var_N(i), Var_N_m(i)] = get_data(a(i),b,n)
    end
    
    properties = gda()
            
    // x labels default
    properties.x_label.text="a";
    properties.x_label.font_style = 6;
    properties.x_label.font_size = 4;
    properties.x_label.foreground = 5;
    properties.x_location = "origin";
    
    // y labels default
    properties.y_label.font_style = 6;
    properties.y_label.font_size = 4;
    properties.y_label.foreground = 3;
    properties.y_location = "origin";
    properties.thickness = 2;
    properties.foreground = 2;
    
    // y name
    properties.y_label.text="$\pi$";
    clf()     
    plot2d(a, Pi, 3)
    xs2eps(gcf(),'Pi.eps')
    
    // y name
    properties.y_label.text="N";
    clf()    
    plot2d(a, N, 3)
    xs2eps(gcf(),'N.eps')
    
    // y name
    properties.y_label.text="Var_N";
    clf()     
    plot2d(a, Var_N, 3)
    xs2eps(gcf(),'Var_N.eps')
    
    // y name
    properties.y_label.text="N*";
    clf()     
    plot2d(a, N_m, 3)
    xs2eps(gcf(),'N_m.eps')
    
    // y name
    properties.y_label.text="Var*_N";
    clf()     
    plot2d(a, Var_N_m, 3)
    xs2eps(gcf(),'Var_N_m.eps')
endfunction

// Сохранение информации в файл *.csv
function[]=save_data(a,b,n)
    [Pi, N, N_m, Var_N, Var_N_m] = get_data(a,b,n)
    data = [a, b, n, Pi, N, Var_N, N_m, Var_N_m]
    write_csv(data, "data.csv", ";")
endfunction

// Получение значений p, p*
function[p, p_m]=get_stats(a,b,n)    
    a_ = 1 - a 
    b_ = 1 - b
    
    // p
    // 1 шаг. p_ij
    for i=1:n+1
        p(i,1) = a_ * b^(i-1)
    end
    
    for i=2:n+1
        for j=2:min(i,n)
            p(i,j) = a_ * get_C(i-1,j-1) * b_^(j-1) * b^(i-j) + ..
                     a * get_C(i-1,j-2) * b_^(j-2) * b^(i-j+1)
        end
    end 
    
    for i=1:n
        p(i,i+1) = a * b_^(i-1)
    end  
    
    p(n+1,n+1) = b_^n + a * get_C(n, 1) * b_^(n-1) * b
    
    // 2 шаг. P_ij
    for i=2:n+1
        for j=1:i-1 
            P(i,j) = sum(p(i,1:j))
        end
    end 
    
    // 3 шаг. q
    q(n+1) = 1
    for i=n:-1:1
        summa = 0
        for j=i+1:n+1
            summa = summa + q(j) * P(j,i)
        end 
        q(i) = summa / p(i,i+1)
    end
    
    // 4 шаг
    clear p
    p(n+1) = 1 / (1 + sum(q(1:n)))
    for i=1:n
        p(i) = p(n+1) * q(i)
    end 
    
    clear q
    clear P
    // p*
    // 1 шаг. p*_ji
    for j=1:n+1
        for i=1:j
            p_mm(j,i) = get_C(j-1, i-1) * b_^(i-1) * b^(j-i)
        end
    end 
    
    // 2 шаг
    for i=1:n+1
        p_m(i) = 0
        for j = i:n+1
            p_m(i) = p_m(i) + p(j)*p_mm(j,i)     
        end
    end
endfunction

// Получение значений Pi, N, Var_N, N*, Var*_N
function[Pi, N, N_m, Var_N, Var_N_m]=get_data(a,b,n)    
    [p,p_m] = get_stats(a,b,n)
        
    // Pi
    Pi = p_m(n+1)
        
    // N
    N=0
    for i=1:n+1
        N = N + (i-1)*p(i)
    end
    
    // N*
    N_m=0
    for i=1:n+1
        N_m = N_m + (i-1)*p_m(i)
    end
    
    // Var_N
    Var_N=0
    for i=1:n+1
        Var_N = Var_N + p(i)*(i-1-N)^2
    end
    
    // Var*_N
    Var_N_m=0
    for i=1:n+1
        Var_N_m = Var_N_m + p_m(i)*(i-1-N_m)^2
    end
endfunction  

// C - число сочетаний без повторений
function[R]=get_C(n,k)
    R = factorial(n) / (factorial(k) * factorial(n-k))
endfunction
