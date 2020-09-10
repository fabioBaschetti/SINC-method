function [con,aon,put] = SINC_fastFT(S,t,K,params_rHeston,j,xi,Xc,Xm,pow)

    % OUTPUT: CoN, AoN and Plain Vanilla put prices from the FFT-version
    % of the SINC formula. This uses N/2 evaluations of the CF for CoN 
    % and AoN options --> N for put options

    if isrow(K)
        % do nothing
    else
        K = K';
    end

    % number of Fourier modes for put prices
    N   = 2^pow;
    
    % construct appropriate truncation range
    Xl  = -Xc+Xm;          
    Xh  =  Xc+Xm;          
    
    % length of the truncation range
    wid =  Xh-Xl;           % = 2*Xc
    
    % sampled Fourier modes \omega_n for CoN
    n   = 1:2:N;      
    wn  = n/wid;
    
    % weights for positive odd frequencies
    Il  = 2 ./ n;
    
    % CF at points wn  for the CoN 
    f1l = exp(-1i*2*pi*wn*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*wn,t);
    % and wn-1i/(2*pi) for the AoN (exponential shift included)
    f2l = exp(-1i*2*pi*wn*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*(wn-1i/(2*pi)),t);
    
    % set even frequencies at zero
    q1l = zeros(1,N);
    q1l(2:2:end) = f1l .* Il;
    q2l = zeros(1,N);
    q2l(2:2:end) = f2l .* Il;
    
    % weights for negative odd frequencies
    Iu  = -flip(Il);
    
    % CF at points  -wn   for the CoN 
    f1u = flip(conj(f1l));
    % and -[wn-1i/(2*pi)] for the AoN (exponential shift included)
    f2u = flip(conj(f2l));
    
    % set even frequencies at zero
    q1u = zeros(1,N);
    q1u(1:2:end) = f1u .* Iu;
    q2u = zeros(1,N);
    q2u(1:2:end) = f2u .* Iu;

    % rearrange positive and negative frequencies...
    q1  = 1i / (2*pi) * [q1l 0 q1u(1:(N-1))];
    q2  = 1i / (2*pi) * [q2l 0 q2u(1:(N-1))];

    % ... for the fft algorithm...
    d1  = fft(q1) + 0.5;
    d2  = fft(q2) + 0.5;
    
    % ... and for computing digital options
    d11 = d1(1:N);
    d12 = d1(N+1:2*N);
    d21 = d2(1:N);
    d22 = d2(N+1:2*N);

    dd1 = [d12 d11];
    dd2 = [d22 d21];

    % define a strike grid...
    m   = (-N):(N-1);
    % ... and compute its buckets (shift Xm included)
    km  = m * wid / (2*N) + Xm;
    Km  = S * exp(km);  
    
    % price CoN and Aon over the buckets
    con = real(Km .* dd1);
    aon = real(S * dd2);
    
    % interpolate on the desired strikes
    con = interp1(Km,con,K,'linear');
    aon = interp1(Km,aon,K,'linear');
    
    % put price as a linear combination of digital options
    put = con - aon;        % eq. (1)
    
    put = put(:);
    
    con = con./K;
    aon = aon /S;
    
    con = con(:);
    aon = aon(:);
    
end