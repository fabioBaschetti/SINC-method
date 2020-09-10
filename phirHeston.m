function out = phirHeston(params_rHeston,xi,n,u,t)
    
    H   = params_rHeston(1);
    rho = params_rHeston(3);
    nu  = params_rHeston(2);

    alpha = H + 0.5;
    
    ti  = (0:(n-1)) / n * t;
    y   = nu * ti.^alpha;
    
    [rm,rp,p1,p2,p3,q1,q2,q3] = dh_Pade33_coeff(alpha,rho,u);
    
    m = size(u,2);
    
    charFun = zeros(1,m);
    for j = 1:m
        dah = dh_Pade33(y,rm(1,j),rp(1,j),p1(1,j),p2(1,j),p3(1,j),q1(1,j),q2(1,j),q3(1,j)); 
        charFun(1,j) = exp(flip(xi) * dah.' * t/n);
    end

    out = charFun;
    
end