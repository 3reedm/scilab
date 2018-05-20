// Сохранение графиков
function[]=save_graphics(b,n,r)
    a = 0.01:0.01:b-0.01
    l = length(a)
    
    for i=1:l
        [Pi(i), N(i), N_m(i), Var_N(i), ..
         Var_N_m(i), w(i), D_w(i)] = get_data(a(i),b,n,r)
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
    
    // y name
    properties.y_label.text="w";
    clf()     
    plot2d(a, w, 3)
    xs2eps(gcf(),'w.eps')
    
    // y name
    properties.y_label.text="D_w";
    clf()     
    plot2d(a, D_w, 3)
    xs2eps(gcf(),'D_w.eps')
endfunction

// Сохранение информации в файл *.csv
function[]=save_data(a,b,n,r)
    [Pi, N, N_m, Var_N, Var_N_m, w, D_w] = get_data(a,b,n,r)
    data = [a, b, n, r, Pi, N, Var_N, N_m, Var_N_m, w, D_w]
    write_csv(data, "data.csv", ";")
endfunction

// Получение значений p, p*, wi
function[p, p_m, wi]=get_stats(a,b,n,r)    
    a_ = 1 - a 
    b_ = 1 - b
    R = n + r
    
    // p
    // 1 шаг. p_ij
    for i=1:n
        p(i,1) = a_ * b^(i-1)
    end
    
    for i=2:n
        for j=2:i
            p(i,j) = a_ * get_C(i-1,j-1) * b_^(j-1) * b^(i-j) + ..
                     a * get_C(i-1,j-2) * b_^(j-2) * b^(i-j+1)
        end
    end 
    
    for i=1:n
        p(i,i+1) = a * b_^(i-1)
    end  
    
    for i=n+1:R+1
        p(i,i-n) = a_ * b^(n+1)
    end  
    
    for i=n+1:R+1
        for j=i-n:min(i,R)
            p(i,j) = a_ * get_C(n+1,i-j) * b_^(n+1-i+j) * b^(i-j) + ..
                     a * get_C(n+1,i-j+1) * b_^(n-i+j) * b^(i-j+1)
        end
    end 
    
    for i=n+1:R
        p(i,i+1) = a * b_^(n+1)
    end  
    
    p(R+1,R+1) = b_^(n+1) + a * n * b_^n * b
    
    // 2 шаг. q
    q(R+1) = 1
    for i=R+1:-1:r+2
        summa = 0
        for j=i+1:R+1
            summa = summa + q(j) * p(j,i)
        end 
        q(i-1) = 1/p(i-1,i) * ((1-p(i,i))*q(i) - summa)
    end
    
    for i=r+1:-1:2
        summa = 0
        for j=i+1:i+n
            summa = summa + q(j) * p(j,i)
        end 
        q(i-1) = 1/p(i-1,i) * ((1-p(i,i))*q(i) - summa)
    end
    
    // 3 шаг
    clear p
    p(R+1) = 1 / (1 + sum(q(1:R)))
    for i=1:R
        p(i) = p(R+1) * q(i)
    end 
    
    clear q
    // p*
    // 1 шаг. p*_ij
    for i=1:n
        for j=1:i
            p_mm(i,j) = get_C(i-1, j-1) * b_^(j-1) * b^(i-j)
        end
    end 
    
    for i=n+1:R+1
        for j=i-n:i
            p_mm(i,j) = get_C(n+1, i-j) * b_^(n+1-i+j) * b^(i-j)
        end
    end 
    
    // 2 шаг
    for i=1:r+1
        p_m(i) = 0
        for j = i:i+n
            p_m(i) = p_m(i) + p(j)*p_mm(j,i)     
        end
    end
    
    for i=r+2:R+1
        p_m(i) = 0
        for j = i:R+1
            p_m(i) = p_m(i) + p(j)*p_mm(j,i)     
        end
    end
    
    clear p_mm
    t = R+1
    // wi
    // 1 шаг. w*_ij
    for j=1:t
        if j<=n+1 then
            w_mm(2,j) = get_C(n, j-1) * b^(j-1) * b_^(n+1-j)
        else
            w_mm(2,j) = 0
        end
    end 
    
    for i=3:t
        for j=1:t
            if j<=(i-1)*n+1 then
                w_mm(i,j) = get_C((i-1)*n, j-1) * b^(j-1) * b_^(n*(i-1)-j+1)
            else
                w_mm(i,j) = 0
            end
        end
    end 
    
    // 2 шаг w_ij
    for j=1:t-1
        w_m(2,j) = 0
        for k=j+1:t
            w_m(2,j) = w_m(2,j) + w_mm(2,k)     
        end
    end
    
    for i=3:t
        for j=1:t
            w_m(i,j) = 0
            for k=1:j-1
                w_m(i,j) = w_m(i,j) + w_mm(i-1,k)*w_m(2,j-k)     
            end
        end
    end
    
    // Pi
    Pi = p_m(R+1)
    // 3 шаг 
    clear w_mm
    wi(1) = 0
    for j=1:n
        wi(1) = wi(1) + p_m(j)     
    end
    wi(1) = 1/(1-Pi) * wi(1)
    
    for i=2:R+1
        wi(i) = 0
        for j=n+1:R
            wi(i) = wi(i) + p_m(j)*w_m(i,j)     
        end
        wi(i) = 1/(1-Pi) * wi(i)
    end
endfunction

// Получение значений Pi, N, Var_N, N*, Var*_N
function[Pi, N, N_m, Var_N, Var_N_m, w, D_w]=get_data(a,b,n,r)    
    R = n + r
    [p,p_m,wi] = get_stats(a,b,n,r)
        
    // Pi
    Pi = p_m(R+1)
        
    // N
    N=0
    for i=1:R+1
        N = N + (i-1)*p(i)
    end
    
    // N*
    N_m=0
    for i=1:R+1
        N_m = N_m + (i-1)*p_m(i)
    end
    
    // Var_N
    Var_N=0
    for i=1:R+1
        Var_N = Var_N + p(i)*(i-1-N)^2
    end
    
    // Var*_N
    Var_N_m=0
    for i=1:R+1
        Var_N_m = Var_N_m + p_m(i)*(i-1-N_m)^2
    end
    
    t = R + 1
    // w
    w=0
    for i=1:t
        w = w + (i-1)*wi(i)
    end
    
    // D_w
    D_w=0
    for i=1:t
        D_w = D_w + wi(i)*(i-1-w)^2
    end
endfunction  

// C - число сочетаний без повторений
function[R]=get_C(n,k)
    R = factorial(n) / (factorial(k) * factorial(n-k))
endfunction
