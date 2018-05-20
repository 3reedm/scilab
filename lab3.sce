// Сохранение графиков
function[]=save_graphics(b,g)
    a = 0.001:0.001:b-0.001
    l = length(a)
    
    for i=1:l
        [Q_1(i), Q_2(i), Q_3(i), Q_0_m(i), Q_1_m(i), w(i)] = get_data(a(i),b,g)
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
    properties.y_label.text="Q_1";
    clf()     
    plot2d(a, Q_1, 3)
    xs2eps(gcf(),'Q_1.eps')
    
    // y name
    properties.y_label.text="Q_2";
    clf()    
    plot2d(a, Q_2, 3)
    xs2eps(gcf(),'Q_2.eps')
    
    // y name
    properties.y_label.text="Q_3";
    clf()     
    plot2d(a, Q_3, 3)
    xs2eps(gcf(),'Q_3.eps')
    
    // y name
    properties.y_label.text="Q_0_m";
    clf()      
    plot2d(a, Q_0_m, 3)
    xs2eps(gcf(),'Q_0_m.eps')
    
    // y name
    properties.y_label.text="Q_1_m";
    clf()     
    plot2d(a, Q_1_m, 3)
    xs2eps(gcf(),'Q_1_m.eps')
    
    // y name
    properties.y_label.text="w";
    clf()     
    plot2d(a, w, 3)
    xs2eps(gcf(),'w.eps')
endfunction

// Сохранение информации в файл *.csv
function[]=save_data(a,b,g)
    [Q_1, Q_2, Q_3, Q_0_m, Q_1_m, w] = get_data(a,b,g)
    data = [a, b, g, Q_1, Q_2, Q_3, Q_0_m, Q_1_m, w]
    write_csv(data, "data.csv", ";")
endfunction

// Получение значений p_0, p_1, p_check_00, p_m, v_0, v_1
function[p_0, p_1, p_check_00, p_m, v_0, v_1]=get_stats(a,b,g)    
    n = 1500
    ro = a / b
    a_ = 1 - a 
    b_ = 1 - b
    g_ = 1 - g
    
    // q_1
    q_0(1) = 1
    q_1(1) = a * q_0(1) / (a_ * b)
    for i=2:n
        up = a * b_ * (1-a_*g_^(i-1))
        down = a_ * b * (1-g_^(i-1))
        q_1(i) = up * q_1(i-1) / down 
    end
    
    // q_0
    for i=2:n
        up = a * b_ * g_^(i-1)
        down = 1 - g_^(i-1)
        q_0(i) = up * q_1(i-1) / down 
    end
    
    // p_0
    p_0(1) = 1 / (sum(q_0) + sum(q_1))
    for i=2:n
        p_0(i) = p_0(1) * q_0(i)
    end  
    
    // p_1
    for i=1:n
        p_1(i) = p_0(1) * q_1(i)
    end 
    
    clear q_0
    clear q_1
    
    // p_check_00
    tmp = 1
    for i=1:n-1
        up = a_ * (b-a*b_*g_^i)
        down = a_*b - a*b_*g_^i
        tmp = tmp * up / down
    end  
    p_check_00 = (1 - ro) / tmp
    
    // p_m
    P_1 = ro
    P_0 = 1 - ro
    p_m(1) = P_0 + P_1*b
    for i=2:n
        p_m(i) = p_1(i-1) * b_
    end
    
    // v_0
    P_1 = ro
    P_0 = 1 - ro
    p_m(1) = P_0 + P_1*b
    for i=2:n
        p_m(i) = p_1(i-1) * b_
    end
    
    // v_1
    for i=1:n
        v_0(2,i) = a_ * (1-g_^i) / i
        v_1(2,i) = v_0(2,i) * b
    end
    
    for j=3:n
        v_0(j,1) = a*v_1(j-1,1) + a_*g_*v_0(j-1,1)
        v_1(j,1) = a*(b_*v_1(j-1,2)+b*v_1(j-1,1)) + ..
                   a_*(b_*v_1(j-1,1)+b*g_*v_0(j-1,1))
        for i=2:n
            v_0(j,i) = a*v_1(j-1,i) + a_*((g_^i)*v_0(j-1,i)+..
                                           (i-1)*(1-g_^i)*v_1(j-1,i-1)/i)
            if (i == n)
                v_1(j,i) = 0
            else
                v_1(j,i) = v_0(j,i)*b + a*b_*v_1(j-1,i+1) + a_*b_*v_1(j-1,i)
            end
        end
    end
endfunction

// Получение значений Q_1, Q_2, Q_3, Q_0_m, Q_1_m, w
function[Q_1, Q_2, Q_3, Q_0_m, Q_1_m, w]=get_data(a,b,g)    
    n = 1500
    ro = a / b
    a_ = 1 - a 
    b_ = 1 - b
    g_ = 1 - g
    
    [p_0, p_1, p_check_00, p_m, v_0, v_1] = get_stats(a,b,g)
    P_0_ = numderivative(list(get_P_0,a,b,a_,b_,g_,p_0(1)), 1)
    P_1_ = numderivative(list(get_P_1,a,b,a_,b_,g_,p_0(1)), 1)    
    
    clear q_0
    clear q_1
    // Q_1
    Q_1 = P_0_ + P_1_
          
    // Q_2
    Q_2 = a*a*b_/(b*(b-a)) + b*P_0_/(b-a)
                                           
    // Q_3
    Q_3 = 0
    for i=2:n
        Q_3 = Q_3 + (i-1) * (p_0(i)+p_1(i))
    end
    
    // Q_0_m
    Q_0_m = (P_0_+b*P_1_) / p_m(1)
    
    // Q_1_m
    Q_1_m = b_ * P_1_ / (1-p_m(1))
    
    // w
    w = 0
    for j=2:n
        for i=1:n
            w = w + j * (v_0(j,i)+v_1(j,i))
        end
    end
endfunction  

// P_0(z)
function[R]=get_P_0(z,a,b,a_,b_,g_,p_0)
    n = 1500
    ro = a / b
    tmp = 1
    for i=1:n-1
       up = a * a * b_ * g_^i * z
       down = a_*b - a*b_*g_^i*z
       tmp = tmp * (1 + up / down)
    end
    R = p_0 * tmp
endfunction

// P_1(z)
function[R]=get_P_1(z,a,b,a_,b_,g_,p_0)
    R = a * get_P_0(z,a,b,a_,b_,g_,p_0) / (a_*b-a*b_*z)
endfunction
