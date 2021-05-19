function out = COSm_discFT(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,N,CP,type)

    if isrow(K)
        % do nothing
    else
        K = K';
    end

    K1  = exp(-(IR-DY)*t)*K;
    
    x   = log(S./K1);
    
    a   = -Xc+Xm;
    b   =  Xc+Xm;

    k   = 0:N-1;
    
    c = a;
    d = 0;
    
    u   = k*pi/(b-a);
    
    phi = charfun(model,params,xi,j,u,t);
    phi(1) = 1;
    m   = length(x);
    eto = zeros(m,N);
    for count = 1:m
        eto(count,:) = exp(1i*pi*k*(x(1,count)-a)/(b-a));
    end
    
    wgt = ones(1,N);
    wgt(1) = 0.5;
    
    psih = d-c;
    psit = (sin(k*pi*(d-a)/(b-a)) - sin(k*pi*(c-a)/(b-a))) .* ((b-a)./(k*pi));
    psi  = [psih psit(2:end)];
    Vk1  = zeros(m,N);
    for count = 1:m
        Vk1(count,:)  = 2/(b-a) * K(1,count) * psi;
    end
    
    conP  = zeros(1,m);
    for count = 1:m
        conP(1,count) = real(wgt.*phi.*eto(count,:))*Vk1(count,:)';
    end
    conP = conP * exp(-IR*t);
    
    add1 = cos(k*pi*(d-a)/(b-a))*exp(d) - cos(k*pi*(c-a)/(b-a))*exp(c); 
    add2 = ((k*pi)/(b-a)) .* (sin(k*pi*(d-a)/(b-a))*exp(d) - sin(k*pi*(c-a)/(b-a))*exp(c));
    fatt = 1./(1+((k*pi)/(b-a)).^2);
    chi  = fatt .* (add1+add2);
    Vk2  = zeros(m,N);
    for count = 1:m
        Vk2(count,:)  = 2/(b-a) * K1(1,count) * chi;
    end
    
    aonP  = zeros(1,m);
    for count = 1:m
        aonP(1,count) = real(wgt.*phi.*eto(count,:))*Vk2(count,:)';
    end
    aonP = aonP * exp(-DY*t);
    
    put = conP - aonP;
    
    call = put+S*exp(-DY*t)-K*exp(-IR*t);
    
    conC = K*exp(-IR*t) - conP;
    aonC = S*exp(-DY*t) - aonP;
    
    if CP*type == 1
        out = call(:);
    elseif CP*type == -1
        out = put(:);
    elseif CP*type == 2
        out = conC(:);
    elseif CP*type == -2
        out = conP(:);
    elseif CP*type == 3
        out = aonC(:);
    elseif CP*type == -3
        out = aonP(:);
    end
    
end