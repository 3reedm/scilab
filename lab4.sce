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
        [N1(i), M1(i), M3(i)] = get_data(a(i),b,z)
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
    
    // y name
    properties.y_label.text="M1";
    clf()    
    plot2d(a, M1, 3)
    xs2eps(gcf(),'M1.eps')
    
    // y name
    properties.y_label.text="M3";
    clf()     
    plot2d(a, M3, 3)
    xs2eps(gcf(),'M3.eps')
endfunction

// Сохранение информации в файл *.csv
function[]=save_data(a,b,z)
    b_ = 0
    for i=2:length(b)
        b_ = b_ + (i-1)*b(i)
    end
    [N1, M1, M3, N1m, N2, N3, N4, N2m, N3m, M2] = get_data(a,b,z)
    data = [a, b, M1, M2, M3, N1, N2, N3, N4, N1m, N2m, N3m]
    write_csv(data, "data.csv", ";")
endfunction

// Получение значений pv, Pzv, bz, pp, Pzp, p, pcheck, Pz, pm, Pzm, 
// Pzmcheck, w, wh, wz
function[pv, Pzv, bz, pp, Pzp, p, pcheck, ..
         Pz, pm, Pzm, Pzmcheck, w, wh, wz]=get_stats(a,b,z)    
    n = 100000
    b_ = 0
    l = length(b)
    for i=2:l
        b_ = b_ + (i-1)*b(i)
    end
    ro = a * b_
    a_ = 1 - a 
    for i=1:l
        B(i) = 0
        for j=i:l
            B(i) = B(i) + b(j)
        end
    end
    
    // ~p
    for i=1:l
        pv(i) = 0
        for j=i:l
            pv(i) = pv(i) + get_C(j-1,i-1)*a^(i-1)*a_^(j-i)*b(j)
        end
    end
    
    // ~P(z)
    Pzv = get_Pzv(z,pv)

    // b(z)
    bz = get_bz(z,b)
    
    // pp
    // Pv
    for i=2:l
        Pv(i) = 0
        for j=i:l
            Pv(i) = Pv(i) + pv(j)
        end
    end
    
    // qp
    for i=2:l
        qp(i) = 0
        for j=2:i-1
            qp(i) = qp(i) + qp(j)*Pv(i-j+2)
        end
        qp(i) = 1/pv(1) * (Pv(i)+qp(i))
    end
    clear Pv
    
    // p+(0)
    pp(1) = 0
    for i=2:l
        pp(1) = pp(1) + qp(i)
    end
    pp(1) = 1 / (1+pp(1))
    
    // p+(i)
    for i=2:l
        pp(i) = pp(1)*qp(i)
    end
    clear qp
    
    // P+(z)
    Pzp = get_Pzp(z,a,b_,pv)
    
    // p
    p(1) = pp(1)
    ph(2) = pp(1) + pp(2)
    for i=3:l
        ph(i) = pp(i)
    end
    
    for i=2:l
        tmp1 = 0
        for k=1:l-1
            tmp2 = 0
            for j=max(2,i-k+1):i
               tmp2 = tmp2 + ph(j)*B(k+1)*get_C(k-1,i-j)*a^(i-j)*a_^(k-i+j-1)
            end
            tmp1 = tmp1 + tmp2
        end
        p(i) = a * tmp1 
    end
    clear ph
    
    // pcheck
    //pcheck(1) = pp(1)
    for i=2:l
        // left part
        lp = 0
        for k=i-1:l-1
            lp = lp + B(k+1)*get_C(k-1, i-2)*a^(i-2)*a_^(k-i+1)
        end
        lp = pp(1) * lp
        
        // right part
        rp = 0
        for j=2:i
            tmp = 0
            for k=i-j+1:l-1
                tmp = tmp + B(k+1)*get_C(k-1, i-j)*a^(i-j)*a_^(k-i+j-1)
            end
            rp = rp + pp(j) * tmp
        end
        
        pcheck(i) = a * (lp+rp)
    end
    
    // P(z)
    Pz = get_Pz(z,a,a_,b,b_,pv,pp(1))
    
    // p*
    pm(1) = 1/a_ * pp(1)
    
    for i=2:l
        // left part
        lp = 0
        for k=i:l-1
            lp = lp + B(k+1)*get_C(k-2, i-2)*a^(i-2)*a_^(k-i)
        end
        lp = pp(1) * lp
        
        // left middle part
        lmp = 0
        for j=2:i
            tmp = 0
            for k=i-j+2:l-1
                tmp = tmp + B(k+1)*get_C(k-2, i-j)*a^(i-j)*a_^(k-i+j-2)
            end
            lmp = lmp + pp(j) * tmp
        end
        
        // right middle part
        rmp = 0
        for k=i+1:l
            rmp = rmp + b(k)*get_C(k-2, i-1)*a^(i-1)*a_^(k-i-1)
        end
        rmp = pp(1) * rmp
        
        // right part
        rp = 0
        for j=2:i+1
            tmp = 0
            for k=i-j+3:l
                tmp = tmp + b(k)*get_C(k-2, i-j+1)*a^(i-j+1)*a_^(k-i+j-3)
            end
            if (i<l)
                rp = rp + pp(j) * tmp
            end
        end
        
        pm(i) = a * (lp+lmp+rmp+rp)
    end
    
    // P*(z)
    Pzm = get_Pzm(z,pm)
    
    // P*check(z)
    Pzmcheck = get_Pzmcheck(z,a,a_,b(1),pp(1),pp(2),b,b_,pv)
    
    // w
    w(1) = 1 - ro
    w(2) = a/a_ * w(1)
    
    for i=2:l-1
        w(i+1) = 0
        for j=1:i-2
            w(i+1) = w(i+1) + w(j+1)*a*b(i-j)
        end
        w(i+1) = (w(i)*(1-a*b(2))-w(1)*a*b(i)-w(i+1)) / a_
    end
    
    // w^
    wh(1) = w(1) + w(2) 
    for i=2:l-1
        wh(i) = w(i+1)
    end
    
    // w(z)
    wz = get_wz(z,a,a_,b,w(1))
endfunction

// Получение значений N1, M1, M3, N1m, N2, N3, N4, N2m, N3m, M2
function[N1, M1, M3, N1m, N2, N3, N4, N2m, N3m, M2,l]=get_data(a,b,z)    
    n = 100000
    b_ = 0
    l = length(b)
    for i=2:l
        b_ = b_ + (i-1)*b(i)
    end
    ro = a * b_
    a_ = 1 - a 
    for i=1:l
        B(i) = 0
        for j=i:l
            B(i) = B(i) + b(j)
        end
    end
    
    [pv, Pzv, bz, pp, Pzp, p, pcheck, ..
     Pz, pm, Pzm, Pzmcheck, w, wh, wz]=get_stats(a,b,z) 
    
    // N1
    N1 = numderivative(list(get_Pzp,a,b_,pv),1)
    
    // N2
    N2 = 0
    for i=1:l
        N2 = N2 + (i-1)*pp(i)
    end
    
    // N3
    N3 = numderivative(list(get_Pz,a,a_,b,b_,pv,pp(1)),1)
    
    // N4
    N4 = 0
    for i=1:l
        N4 = N4 + (i-1)*p(i)
    end
    
    // N1m
    N1m = numderivative(list(get_Pzm,pm),1)
    
    // N2m
    N2m = numderivative(list(get_Pzmcheck,a,a_,b(1),pp(1),pp(2),b,b_,pv),1)
    
    // N3m
    N3m = 0
    for i=1:l-1
        N3m = N3m + (i-1)*pm(i)
    end
    
    // M1
    M1 = 0
    for i=1:l
        M1 = M1 + (i-1)*w(i)
    end
    
    // M2
    M2 = numderivative(list(get_wz,a,a_,b,w(1)),1)
    
    // M3
    M3 = 0
    for i=1:l-1
        M3 = M3 + (i-1)*wh(i)
    end
    l= a*(M3+b_)
endfunction  

// C - число сочетаний без повторений
function[R]=get_C(n,k)
    R = factorial(n) / (factorial(k) * factorial(n-k))
endfunction

// ~P(z)
function[R]=get_Pzv(z,pv)
    l = length(pv)
    R = 0
    for i=1:l
        R = R + z^(i-1)*pv(i)
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

// P+(z)
function[R]=get_Pzp(z,a,b_,pv)
    up = (1-z) * get_Pzv(z,pv) * (1-a*b_)
    down = get_Pzv(z,pv) - z
    R = up/down
endfunction

// P(z)
function[R]=get_Pz(z,a,a_,b,b_,pv,pp0)
    bbeta = get_bz(a*z+a_,b)
    Pzp = get_Pzp(z,a,b_,pv)
    up = (1-bbeta) * (Pzp - pp0 + z*pp0)
    down = 1 - z
    R = up/down + pp0
endfunction

// P*(z)
function[R]=get_Pzm(z,pm)
    l = length(pm)
    R = 0
    for i=1:l
        R = R + z^(i-1)*pm(i)
    end
endfunction

// P*check(z)
function[R]=get_Pzmcheck(z,a,a_,b0,pp0,pp1,b,b_,pv)
    bbeta1 = get_bz(a*z+a_,b)
    bbeta2 = get_bz(a_,b)
    Pzp = get_Pzp(z,a,b_,pv)
    tmp = a*z + a_
    up1 = tmp - bbeta1
    down1 = tmp * (1-tmp)
    tmp2 = bbeta1/tmp - bbeta2/a_
    R = pp0/a_ + a*pp0*z*up1/down1 + a*(Pzp-pp0)*up1/down1 + ..
        a*pp0*tmp2 + a/z*(Pzp-pp0)*tmp2 + a/a_*(Pzp-pp0-pp1*z)*(bbeta2-b0)
endfunction

// w(z)
function[R]=get_wz(z,a,a_,b,w0)
    bbeta = get_bz(z,b)
    up = (z-1) * (a_+a*bbeta)
    down = z - a_ - a*bbeta
    R = up*w0/down
endfunction
