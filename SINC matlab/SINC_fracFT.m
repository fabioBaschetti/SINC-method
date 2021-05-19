function out = SINC_fracFT(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,pow,e,CP)

    % OUTPUT: Plain Vanilla put/call prices from the frFFT-version of the 
    % SINC formula. 
    % This uses N/2 evaluations of the CF for CoN and AoN options --> N for put options

    if isrow(K)
        % do nothing
    else
        K = K';
    end

    % number of Fourier modes for Eur prices
    N   = 2^pow;
    
    % construct appropriate truncation range
    Xl  = -Xc+Xm;          
    Xh  =  Xc+Xm;          
    
    % length of the truncation range
    wid =  Xh-Xl;           % = 2*Xc
    
    % sampled Fourier modes \omega_n for CoN put
    n   = 1:2:(N/2);      
    wn  = n/wid;
    
    % weights for positive odd frequencies
    Il  = 2 ./ n;
    
    % CF at points wn  for the CoN put
    f1l = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn   ,t);
    % and wn-1i/(2*pi) for the AoN put (exponential shift included)
    f2l = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn-1i,t);
    
    % set even frequencies at zero
    q1l = zeros(1,N/2);
    q1l(2:2:end) = f1l .* Il;
    q2l = zeros(1,N/2);
    q2l(2:2:end) = f2l .* Il;
    
    % weights for negative odd frequencies
    Iu  = -flip(Il);
    
    % CF at points  -wn   for the CoN put
    f1u = flip(conj(f1l));
    % and -[wn-1i/(2*pi)] for the AoN put (exponential shift included)
    f2u = flip(conj(f2l));
    
    % set even frequencies at zero
    q1u = zeros(1,N/2);
    q1u(1:2:end) = f1u .* Iu;
    q2u = zeros(1,N/2);
    q2u(1:2:end) = f2u .* Iu;

    % rearrange positive and negative frequencies ...
    q1  = 1i / (2*pi) * [0 q1u(1:(N/2-1)) q1l];
    q2  = 1i / (2*pi) * [0 q2u(1:(N/2-1)) q2l];

    % ... for the frFFT algorithm ...
    d1  = frFFT_circ(q1,e) + 0.5;
    d2  = frFFT_circ(q2,e) + 0.5;
    
    dd1 = d1;
    dd2 = d2;

    % define a strike grid...
    m   = (-N/2):(N/2-1);
    % ... and compute its buckets (shift Xm included)
    beta = e / N;
    km  = beta * m * wid + (IR-DY)*t + Xm;
    Km  = S * exp(km);  
    
    % price CoN and Aon put over the buckets
    conP = exp(-IR*t).*real(Km .* dd1);
    aonP = exp(-DY*t).*real(S * dd2);
    
    % interpolate on the desired strikes
    conP = interp1(Km,conP,K,'linear');
    aonP = interp1(Km,aonP,K,'linear');
    
    % put price as a linear combination of digital options
    put = conP - aonP;        % eq. (1)
    
    % call price from put-call parity
    call = put+S*exp(-DY*t)-K*exp(-IR*t);
    
    if CP == 1
        out = call(:);
    elseif CP == -1
        out = put(:);
    end
   
end