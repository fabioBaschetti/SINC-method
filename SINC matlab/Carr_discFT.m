function out = Carr_discFT(S,t,K,IR,DY,model,params,xi,j,dump,ub,N,CP)
    
    if isrow(K)
        % do nothing
    else
        K = K';
    end
    
    eta   = ub/(N-1);
    v     = eta * (0:N-1);
    
    K1    = exp(-(IR-DY)*t)*K;
    k     = log(K1/S);
    
    phi   = charfun(model,params,xi,j,v-(dump+1)*1i,t); 
    psi   = phi ./ (dump^2+dump-v.^2+1i*v*(2*dump+1));
    
    m     = length(k);
    eto   = zeros(m,N);
    for count = 1:m
        eto(count,:) = exp(-1i*v*k(1,count));
    end
    
    % simpson
    w     = 1/3*[1 (3+(-1).^(2:N-1)) 1];
%     % naive
%     w = [0.5 ones(1,N-1)];
    
    tbsum = zeros(m,N);
    for count = 1:m
        tbsum(count,:) = eta*eto(count,:).*psi.*w;
    end
    
    call  = zeros(1,m); 
    for count = 1:m
        call(1,count)  = S*exp(-DY*t)*exp(-dump*k(1,count))/pi * real(sum(tbsum(count,:)));
    end
    
    put = call - S*exp(-DY*t) + K*exp(-IR*t);
    
    if CP == 1
        out = call(:);
    elseif CP == -1
        out = put(:);
    end
    
end