function price = Lewis_quadra(S,t,K,IR,DY,model,params,xi,j,CP)
    price = arrayfun(@(K) LewisCmp(S,t,K,IR,DY,model,params,xi,j,CP),K);
    price = price(:);
end

function out = LewisCmp(S,t,K,IR,DY,model,params,xi,j,CP)
    
    K1  = exp(-(IR-DY)*t)*K;
    k   = log(K1/S);
    integrand = @(u) real(exp(-1i*u*k).*charfun(model,params,xi,j,u-1i/2,t)) / (u.^2+1/4);
    int = integral(integrand,0,Inf,'ArrayValued',true,'RelTol',1e-14);
    put = exp(-DY*t)*(K1 - sqrt(S*K1)/pi*int);
    
    call = put+S*exp(-DY*t) - K*exp(-IR*t);
    
    if CP == 1
        out = call(:);
    elseif CP == -1
        out = put(:);
    end

end