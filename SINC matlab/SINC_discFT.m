function out = SINC_discFT(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,N,CP,type)

    % OUTPUT: CoN, AoN and Plain Vanilla put/call prices from the 
    % FT-version of the SINC formula. 
    % This uses N/2 evaluations of the CF for CoN and AoN options --> N for put/call options
    
    if isrow(K)
        % do nothing
    else
        K = K';
    end
    
    % construct appropriate truncation range
    Xl  = -Xc+Xm;
    Xh  =  Xc+Xm;
    
    % length of the truncation range
    wid =  Xh-Xl;           % = 2*Xc
    
    % sampled Fourier modes \omega_n for CoN put
    n   = 1:2:(N/2);            
    wn  = n/wid;            
    
    % log-strike (shift Xm included)
    k   = log(K/S) - (IR-DY)*t - Xm;
    
    % compute sin and cos in formula (10)
    sn  = zeros(size(k,2),size(wn,2));
    cs  = zeros(size(k,2),size(wn,2));
    
    for i=1:size(k,2)
        sn(i,:) = sin(2*pi*k(i)*wn);
        cs(i,:) = cos(2*pi*k(i)*wn);
    end
    
    % CF at points  wn for the CoN put
    f1  = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn   ,t); 
    % and wn-1i/(2*pi) for the AoN put (exponential shift included)
    f2  = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn-1i,t);
    
    % compute sums in eq. (8) and (9) for CoN and AoN put, resp.ly 
    ad1 = zeros(1,size(k,2));
    ad2 = zeros(1,size(k,2));
    
    for i=1:size(k,2)
        ad1(i) = sum((sn(i,:).*real(f1) - cs(i,:).*imag(f1)) ./ n);
        ad2(i) = sum((sn(i,:).*real(f2) - cs(i,:).*imag(f2)) ./ n);
    end
    
    % conclude CoN and AoN put prices (eq. (8) and (9))
    conP = K*exp(-IR*t).*(0.5+(2/pi)*ad1);
    aonP = S*exp(-DY*t).*(0.5+(2/pi)*ad2);
    
    % put price as a linear combination of digital options
    put = conP - aonP;   % eq. (1)
    
    % call price from put-call parity ...
    call = put+S*exp(-DY*t)-K*exp(-IR*t);
    
    % ... and the same for digital call options
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