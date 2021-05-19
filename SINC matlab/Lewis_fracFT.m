function out = Lewis_fracFT(S,t,K,IR,DY,model,params,xi,j,ub,pow,e,CP)
    
    if isrow(K)
        % do nothing
    else
        K = K';
    end
    
    N     = 2^pow;
    
    eta   = ub/(N-1);
    u     = eta * (0:N-1);
    
    beta  = e / N;
    lambda = 2*pi*beta /eta;
    
    b      = N * lambda / 2;
    
    ku    = -b + lambda*(0:N-1) + log(S) + (IR-DY)*t;
    Ku    = exp(ku);
    
    phi   = charfun(model,params,xi,j,u-0.5*1i,t);
    psi   = phi ./ (u.^2+0.25);
    
%     % simpson
%     w     = 1/3*[1 (3+(-1).^(2:N-1)) 1];
    % naive
    w = [0.5 ones(1,N-1)];
    
    h     = exp(1i*b*u) .* psi * eta .* w;
    x     = real(frFFT_chou(h,e));
    C     = S*exp(-DY*t) - sqrt(S*Ku)/pi * exp(-0.5*(IR+DY)*t) .* x;
        
    Kgrid = Ku';
    [~,I] = sort(abs(bsxfun(@minus,Kgrid,K)));
    Cgrid =  C';
    dC    = Cgrid(I(2,:),:)-Cgrid(I(1,:),:);
    dim   = size(Cgrid,2);
    dK    = repmat(Kgrid(I(2,:))-Kgrid(I(1,:)),1,dim);
    call  = Cgrid(I(1,:),:) + dC./dK.*repmat((K'-Kgrid(I(1,:))),1,dim);
        
    K   = K';
    put = call - S*exp(-DY*t) + K*exp(-IR*t);
    
    if CP == 1
        out = call(:);
    elseif CP == -1
        out = put(:);
    end
    
end