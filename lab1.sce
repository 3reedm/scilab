// Сохранение графиков
function[]=save_graphics(b)
    a = 0.001:0.001:b-0.001
    l = length(a)
    
    for i=1:l
        [N(i), Q(i), N_m(i), Q_m(i), Var_N(i), ..
         Var_N_m(i), w(i), v(i), D_w(i), D_v(i)] = get_data(a(i),b)
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
    properties.y_label.text="N";
    clf()     
    plot2d(a, N, 3)
    xs2eps(gcf(),'N.eps')
    
    // y name
    properties.y_label.text="Q";
    clf()    
    plot2d(a, Q, 3)
    xs2eps(gcf(),'Q.eps')
    
    // y name
    properties.y_label.text="N*";
    clf()     
    plot2d(a, N_m, 3)
    xs2eps(gcf(),'N_m.eps')
    
    // y name
    properties.y_label.text="Q*";
    clf()      
    plot2d(a, Q_m, 3)
    xs2eps(gcf(),'Q_m.eps')
    
    // y name
    properties.y_label.text="Var_N";
    clf()     
    plot2d(a, Var_N, 3)
    xs2eps(gcf(),'Var_N.eps')
    
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
    properties.y_label.text="v";
    clf()     
    plot2d(a, v, 3)
    xs2eps(gcf(),'v.eps')
    
    // y name
    properties.y_label.text="D_w";
    clf()      
    plot2d(a, D_w, 3)
    xs2eps(gcf(),'D_w.eps')
    
    // y name
    properties.y_label.text="D_v";
    clf()     
    plot2d(a, D_v, 3)
    xs2eps(gcf(),'D_v.eps')
endfunction

// Сохранение информации в файл *.csv
function[]=save_data(a,b)
    [N, Q, N_m, Q_m, Var_N, Var_N_m, w, v, D_w, D_v] = get_data(a,b)
    data = [a, b, N, Var_N, Q, N_m, ..
            Var_N_m, Q_m, w, v, D_w, D_v]
    write_csv(data, "data.csv", ";")
endfunction

// Получение значений p, p*, v_k, w_k
function[p, p_m, v_k, w_k]=get_stats(a,b)    
    n = 100
    k = 100
    ro = a/b
    a_ = 1 - a 
    b_ = 1 - b
    delta = a*b_ / (a_*b)
    
    // p
    p(1) = 1 - ro
    for i=2:n
        p(i) = 1/b_ * delta^(i-1) * p(1)
    end
    
    // p*
    p_m(1) = 1 - delta
    for i=2:n
        p_m(i) = delta^(i-1) * p_m(1)
    end  
    
    // v_k
    v_k(1) = 0
    for i=2:k
        v_k(i) = (1-b_/a_) * (b_/a_)^(i-2)
    end  
    
    // w_k
    w_k(1) = p_m(1)
    for i=2:k
        w_k(i) = (1-p_m(1)) * (1-b_/a_) * (b_/a_)^(i-2)
    end
endfunction

// Получение значений N, Q, N*, Q*, Var_N, Var*_N, w, v, D_w, D_v
function[N, Q, N_m, Q_m, Var_N, Var_N_m, w, v, D_w, D_v]=get_data(a,b)    
    ro = a/b
    a_ = 1 - a 
    b_ = 1 - b
    delta = a*b_ / (a_*b)
        
    // N
    N = a_*ro / (1-ro)
    
    // Q
    Q = b_*ro^2 / (1-ro)
    
    // N*
    N_m = b_*ro / (1-ro)
    
    // Q*
    Q_m = (b_*ro)^2 / (a_*(1-ro))
    
    // p(1)
    p = 1 - ro
    // Var_N
    [J,H] = numderivative(list(get_P,a,b,a_,b_,p), 1)
    Var_N = H + N - N^2
    
    // Var*_N
    [J,H] = numderivative(list(get_P_m,delta), 1)
    Var_N_m = H + N_m - N_m^2
    
    // w
    w = a*b_ / (b*(a_-b_))
    
    // v
    v = w + 1/b
    
    // p*(1)
    p_m = 1 - delta
    // D_w
    [J,H] = numderivative(list(get_o,a_,b_,p_m), 1)
    D_w = H + w - w^2
    
    // D_v
    [J,H] = numderivative(list(get_fi,a_,b_), 1)
    D_v = H + v - v^2
endfunction  

// P(z)
function[R]=get_P(z,a,b,a_,b_,p0)
    up = a_*b + z*a*b
    down = a_*b - z*a*b_
    R = up/down * p0
endfunction

// P*(z)
function[R]=get_P_m(z,delta)
    R = (1-delta) / (1-z*delta)
endfunction

// o(z)
function[R]=get_o(z,a_,b_,p_m0)
    up = (a_-b_) * z
    down = a_ - b_*z
    R = p_m0 + (1-p_m0)*up/down
endfunction

// fi(z)
function[R]=get_fi(z,a_,b_)
    up = (a_-b_) * z
    down = a_ - b_*z
    R = up/down
endfunction
