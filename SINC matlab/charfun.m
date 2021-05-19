function phi = charfun(model,params,xi,n,u,t)

    if strcmp(model,'GBM')
        
        sigma = params(1);
        phi = exp(-0.5*1i*u*t*sigma^2 - 0.5*t*(u*sigma).^2);
    
    elseif strcmp(model,'Heston')
        
        lambda = params(1);
        v_bar  = params(2);
        eta    = params(3);
        v0     = params(4);
        rho    = params(5);
        
        beta  = lambda - rho*eta*1i*u;
        alpha = -0.5*u.^2 - 0.5*1i*u;
        zeta  = 0.5*eta^2;

        d = sqrt(beta.^2-4*alpha*zeta);

        rm = (beta-d)/eta^2;
        rp = (beta+d)/eta^2;

        g = rm./rp;

        C = rm*t - 2/eta^2*log((1-g.*exp(-d*t))./(1-g));
        D = rm.*(1-exp(-d*t))./(1-g.*exp(-d*t));

        phi = exp(v_bar*lambda*C + v0*D);
        
    elseif strcmp(model,'VG')
        
        sigma = params(1);
        theta = params(2);
        nu    = params(3);      
        
        fac1 = exp((1i*u*t)/nu * log(1-nu*theta-0.5*(nu*sigma^2)));
        fac2 = (1./(1+2*nu*(0.5*u*sigma).^2-1i*u*nu*theta)).^(t/nu);
        
        phi = fac1 .* fac2;
        
    elseif strcmp(model,'CGMY')
    
        C = params(1);
        G = params(2);
        M = params(3);
        Y = params(4);
        
        g = @(u) gamma(-Y)*((M-1i*u).^Y - M^Y + (G+1i*u).^Y - G^Y);
        phi = exp(C*t*(g(u)-1i*u*g(-1i)));
        
    elseif strcmp(model,'rough Heston')
        
        H   = params(1);
        rho = params(3);
        nu  = params(2);

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

        phi = charFun;
        
    end

end
