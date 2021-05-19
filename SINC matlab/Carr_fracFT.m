function out = Carr_fracFT(S,t,K,IR,DY,model,params,xi,j,dump,ub,pow,e,CP)

    if isrow(K)
        % do nothing
    else
        K = K';
    end

    N      = 2^pow;
    
    x0     = log(S);
    
    eta    = ub/(N-1);
    v      = eta *(0:N-1);
    
    beta   = e / N;
    lambda = 2*pi*beta /eta;
    
    b      = N * lambda/2;
    
    ku     = lambda*(0:N-1) - b + x0;

    u      = v-1i*(dump+1);
    phi    = charfun(model,params,xi,j,u,t);
    psi    = exp(-IR*t) * exp(1i*u*(x0+(IR-DY)*t)) .*phi ./ (dump^2+dump-v.^2+1i*(2*dump+1)*v);
    
    % simpson
    w     = 1/3*[1 (3+(-1).^(2:N-1)) 1];
%     % naive
%     w = [0.5 ones(1,N-1)];
    
    h      = exp(v*1i*(b-x0)) .*psi * eta .* w;
    x      = real(frFFT_chou(h,e));
    C      = exp(-dump*ku)/pi .* x;
    
    Ku     = exp(ku);
    
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