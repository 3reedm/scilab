// Сохранение графиков
function[]=save_graphics(b)
    b_ = 0
    for i=2:length(b)
        b_ = b_ + (i-1)*b(i)
    end
    b_ = 1/b_
    a = 0.001:0.001:b_-0.001
    l = length(a)
    
    step = 0
    for i=1:l
        a_ = 1 - a(i)
        tmp = a_ / (1+a_)
        if (get_bz(a_,b) > tmp) then
            step = step + 1
            anew(step) = a(i)
        end
    end
    
    lnew = step
    for i=1:lnew
        [N(i), v1(i), g1_(i)] = get_data(anew(i),b)
    end
    
    properties = gda()
            
    // x labels default
    properties.x_label.text="a*";
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
    plot2d(anew, N, 0)
    xs2eps(gcf(),'N.eps')
    
    // y name
    properties.y_label.text="v1";
    clf()     
    plot2d(anew, v1, 0)
    xs2eps(gcf(),'v1.eps')
    
    // y name
    properties.y_label.text="g1_";
    clf()     
    plot2d(anew, g1_, 0)
    xs2eps(gcf(),'g1_.eps')
endfunction

// Сохранение информации в файл *.csv
function[]=save_data(a,b)
    b_ = 0
    for i=2:length(b)
        b_ = b_ + (i-1)*b(i)
    end
    [N, v1, g1_, s1_, vhm, s2_, g2_, v2] = get_data(a,b)
    data = [a, b, s1_, s2_, N, g1_, g2_, v1, v2, vhm(1,:), vhm(2,:)]
    write_csv(data, "data.csv", ";")
endfunction

// Получение значений p, p(i,j), g(k), v(k), v(k,m)
function[p, pk, g, v, vm]=get_stats(a,b)    
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
    
    // p0
    bbeta = get_bz(a_,b)
    p(1) = 1 - a_*(1-bbeta)/bbeta
    p(2) = 1/bbeta * (1-bbeta) * p(1)
    
    // pk
    pk(2,2) = a/a_ * p(1)
    for i=3:l
        pk(2,i) = a/a_ * (p(i-1)-pk(2,i-1))
        p(i) = 1/bbeta * (1-bbeta) * (p(i-1)-pk(2,i-1))
    end
    
    for i=3:l
        pk(i,2) = 0
        for j=i:l
            pk(i,2) = pk(i,2) + a_^(j-i)*b(j) 
        end
        pk(i,2) = a * (p(1)+p(2)) * pk(i,2) 
    end
    
    for i=3:l
        for k=3:l
            pk(i,k) = 0
            for j=i:l
                pk(i,k) = pk(i,k) + a_^(j-i)*b(j) 
            end
            pk(i,k) = (a*(p(k-1)+pk(2,k-1))+a*p(k)) * pk(i,k)
        end
    end
    
    // g
    g(1) = 0
    for k=2:l
        tmp1 = 0
        tmp2 = 0
        for i=2:k-1
            tmp1 = tmp1 + a_^(i-2)*a*b(i)*g(k-i+1)
            
            tmp3 = 0
            for j=1:k-i
                tmp3 = tmp3 + g(j)*g(k-i-j+2)
            end
            tmp2 = tmp2 + a_^(i-2)*a*B(i+1)*tmp3
        end
        g(k) = a_^(k-1)*b(k) + tmp1 + tmp2
    end
    
    // v
    v(1) = 0
    for k=2:l
        tmp1 = 0
        for i=2:k-1            
            tmp2 = 0
            for j=1:k-i
                tmp2 = tmp2 + g(j)*v(k-i-j+2)
            end
            tmp1 = tmp1 + a_^(i-2)*a*B(i+1)*tmp2
        end
        v(k) = a_^(k-2)*b(k) + tmp1
    end
    
    // vm
    for m=2:l
        vm(1,m) = 0
    end
    
    for m=2:l
        for k=2:m-1
            tmp1 = 0
            for i=2:k-1            
                tmp2 = 0
                for j=1:k-i
                    tmp2 = tmp2 + g(j)*v(k-i-j+2)
                end
                tmp1 = tmp1 + a_^(i-2)*a*tmp2
            end
            vm(k,m) = tmp1
        end
        for k=m:l
            tmp1 = 0
            for i=2:m-1            
                tmp2 = 0
                for j=1:k-i
                    tmp2 = tmp2 + g(j)*v(k-i-j+2)
                end
                tmp1 = tmp1 + a_^(i-2)*a*tmp2
            end
            vm(k,m) = a_^(m-2) + tmp1
        end
    end
endfunction

// Получение значений N, v1, g1_, s1_, vhm
function[N, v1, g1_, s1_, vhm, s2_, g2_, v2]=get_data(a,b)    
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
    
    [p, pk, g, v, vm]=get_stats(a,b)
    
    // s1_
    s1_ = 0
    for i=1:l
        s1_ = s1_ + (i-1)*p(i)
    end
    
    // s2_
    bbeta = get_bz(a_,b)
    s2_ = a_*(1-bbeta) / (a*bbeta)
    
    // N
    N = 0
    for i=2:l
        for j=2:l
            N = N + (i-1)*pk(i,j)
        end
    end
    
    // g1_
    g1_ = numderivative(list(get_yz,a,a_,b),1)
    
    // g2_
    g2_ = 0
    for i=1:l
        g2_ = g2_ + (i-1)*g(i)
    end
    
    // v1
    v1 = numderivative(list(get_vz,a,a_,b),1)
    
    // g2_
    v2 = 0
    for i=1:l
        v2 = v2 + (i-1)*v(i)
    end
    
    // vhm
    for m=2:l
        vhm(1,m) = numderivative(list(get_vzm,a,a_,b,m-1),1)
        
        vhm(2,m) = 0 
        for i=2:l
            vhm(2,m) = vhm(2,m) + (i-1)*vm(i,m)
        end
    end
endfunction  

// b(z)
function[R]=get_bz(z,b)
    l = length(b)
    R = 0
    for i=1:l
        R = R + z^(i-1)*b(i)
    end
endfunction

// y(z)
function[R]=get_yz(z,a,a_,b)
    arg = a_*z
    bbeta = get_bz(arg,b)
    up = (1-arg)*(a_-a*bbeta) - ..
         sqrt((1-arg)^2*(a_-a*bbeta)^2-4*a*(arg-bbeta)*a_*(1-arg)*bbeta)
    down = 2 * a * (arg-bbeta)
    R = up/down
endfunction

// v(z)
function[R]=get_vz(z,a,a_,b)
    arg = a_*z
    bbeta = get_bz(arg,b)
    up = (1-arg) * bbeta
    down = a_*(1-arg) - a*(arg-bbeta)*get_yz(z,a,a_,b)
    R = up/down
endfunction

// v(z;m)
function[R]=get_vzm(z,a,a_,b,m)
    ggamma = get_yz(z,a,a_,b)
    vetta = get_vz(z,a,a_,b)
    arg = a_*z
    R = arg^m/a_ + a*z/a_*(1-arg^(m-1))/(1-arg)*ggamma*vetta
endfunction
