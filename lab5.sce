// Сохранение графиков
function[]=save_graphics(b,z)
    b_ = 0
    for i=2:length(b)
        b_ = b_ + (i-1)*b(i)
    end
    b_ = 1/b_
    a = 0.001:0.001:b_-0.001
    l = length(a)
    
    for i=1:l
        N1(i) = get_data(a(i),b,z)
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
    properties.y_label.text="N1";
    clf()     
    plot2d(a, N1, 3)
    xs2eps(gcf(),'N1.eps')
endfunction

// Сохранение информации в файл *.csv
function[]=save_data(a,b,z)
    b_ = 0
    for i=2:length(b)
        b_ = b_ + (i-1)*b(i)
    end
    [N1, N2, N3, D] = get_data(a,b,z)
    data = [a, b, N1, N2, N3, D]
    write_csv(data, "data.csv", ";")
endfunction

// Получение значений p_n_(m), P(z)
function[p, Pz]=get_stats(a,b,z)    
    n = 10000
    b_ = 0
    l = length(b)
    for i=2:l
        b_ = b_ + (i-1)*b(i)
    end
    ro = a / b_
    a_ = 1 - a 
    for i=1:l
        B(i) = 0
        for j=i:l
            B(i) = B(i) + b(j)
        end
    end
    
    // p_n_(m)
    for i=1:n
        if (i==1) then
            p(i,1) = 1 - a*B(2)
        elseif (i==2) then
            p(i,1) = a * B(2)
        else
            p(i,1) = 0
        end
    end
    
    for i=2:l-1
        p(1,i) = (1-a*B(i+1)) * p(1,i-1)
    end
    
    for i=2:n
        for j=2:l-1
            p(i,j) = (1-a*B(j+1))*p(i,j-1) + a*B(j+1)*p(i-1, j-1)
        end
    end
    
    // P(z)
    Pz = get_Pz(z,a,B)
endfunction

// Получение значений N1, M1, M3, N1m, N2, N3, N4, N2m, N3m, M2
function[N1, N2, N3, D]=get_data(a,b,z)    
    n = 10000
    b_ = 0
    l = length(b)
    for i=2:l
        b_ = b_ + (i-1)*b(i)
    end
    ro = a / b_
    a_ = 1 - a 
    for i=1:l
        B(i) = 0
        for j=i:l
            B(i) = B(i) + b(j)
        end
    end
    
    [p, Pz]=get_stats(a,b,z) 
    
    // N1
    N1 = numderivative(list(get_Pz,a,B),1)
    
    // N2
    N2 = 0
    for i=1:l-1
        N2 = N2 + B(i+1)
    end
    N2 = a * N2
    
    // N3
    N3 = 0
    for i=2:n
        tmp = 0
        for j=3:l-1
            tmp = tmp + p(i,j)
        end
        N3 = N3 + (i-1)*tmp
    end
    
    // D
    [J,H] = numderivative(list(get_Pz,a,B),1)
    D = H + J - J^2
endfunction  

// Получение значений N1, M1, M3, N1m, N2, N3, N4, N2m, N3m, M2
function[N1x, N2x]=get_data_x(a,b,l)    
    n = 10000
    b_ = 0
    k = length(b)
    for i=2:k
        b_ = b_ + (i-1)*b(i)
    end
    ro = a / b_
    a_ = 1 - a 
    for i=1:k
        B(i) = 0
        for j=i:k
            B(i) = B(i) + b(j)
        end
    end
    
    l_ = 0
    kk = length(l)
    for i=2:kk
        l_ = l_ + (i-1)*l(i)
    end
    
    // N1
    N1x = numderivative(list(get_Pzx,l,B),1)    
    // N2
    N2x = 0
    for i=1:k-1
        N2x = N2x + B(i+1)
    end
    N2x = l_ * N2x 
    N1x = N2x
endfunction

// P(z)
function[R]=get_Pz(z,a,B)
    l = length(B)
    R = 1
    for i=1:l-1
        R = R * (1-a*B(i+1)*(1-z))
    end
endfunction

// L(z)
function[R]=get_Lz(z,l)
    k = length(l)
    R = 0
    for i=2:k
        R = R + z^(i-1)*l(i)
    end
endfunction

// P[x](z)
function[R]=get_Pzx(z,l,B)
    k = length(B)
    R = 1
    for i=1:k-1
        R = R * get_Lz(1-(1-z)*B(i+1),l)
    end
endfunction
