function out = Lewis_discFT(S,t,K,IR,DY,model,params,xi,j,ub,N,CP)
    
    if isrow(K)
        % do nothing
    else
        K = K';
    end
    
    eta   = ub/(N-1);
    u     = eta * (0:N-1);
    
    K1    = exp(-(IR-DY)*t)*K;
    k     = log(K1/S);
    
    phi   = charfun(model,params,xi,j,u-0.5*1i,t);
    psi   = phi ./ (u.^2+0.25); 
    
    m     = length(k);
    eto   = zeros(m,N);
    for count = 1:m
        eto(count,:) = exp(-1i*u*k(1,count));
    end
    
%     % simpson
%     w     = 1/3*[1 (3+(-1).^(2:N-1)) 1];
    % naive
    w = [0.5 ones(1,N-1)];
  
    tbsum = zeros(m,N);
    for count = 1:m
        tbsum(count,:) = eta*eto(count,:).*psi.*w;
    end
    
    put  = zeros(1,m); 
    for count = 1:m
        put(1,count)  = exp(-DY*t)*(K1(count) - sqrt(S*K1(count))/pi * real(sum(tbsum(count,:))));
    end
    
    call = put+S*exp(-DY*t) - K*exp(-IR*t);
    
    if CP == 1
        out = call(:);
    elseif CP == -1
        out = put(:);
    end
    
end