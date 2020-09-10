function [con,aon,put] = SINC_discFT(S,t,K,params_rHeston,j,xi,Xc,Xm,N)

    % OUTPUT: CoN, AoN and Plain Vanilla put prices from the FT-version
    % of the SINC formula. This uses N/2 evaluations of the CF for CoN 
    % and AoN options --> N for put options
    
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
    
    % sampled Fourier modes \omega_n for CoN
    n   = 1:2:N;            
    wn  = n/wid;            
    
    % log-strike (shift Xm included)
    k   = log(K/S) - Xm;
    
    % compute sin and cos in formula (10)
    sn  = zeros(size(k,2),size(wn,2));
    cs  = zeros(size(k,2),size(wn,2));
    
    for i=1:size(k,2)
        sn(i,:) = sin(2*pi*k(i)*wn);
        cs(i,:) = cos(2*pi*k(i)*wn);
    end
    
    % CF at points wn  for the CoN 
    f1  = exp(-1i*2*pi*wn*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*wn,t);                 
    % and wn-1i/(2*pi) for the AoN (exponential shift included)
    f2  = exp(-1i*2*pi*wn*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*(wn-1i/(2*pi)),t);
    
    % compute sums in eq. (8) and (9) for CoN and Aon, resp.ly 
    ad1 = zeros(1,size(k,2));
    ad2 = zeros(1,size(k,2));
    
    for i=1:size(k,2)
        ad1(i) = sum((sn(i,:).*real(f1) - cs(i,:).*imag(f1)) ./ n);
        ad2(i) = sum((sn(i,:).*real(f2) - cs(i,:).*imag(f2)) ./ n);
    end
    
    % conclude CoN and AoN put prices (eq. (8) and (9))
    con = (0.5+2/pi*ad1);
    aon = (0.5+2/pi*ad2);
    
    % put price as a linear combination of digital options
    put = K.*con - S*aon;   % eq. (1)
    
    put = put(:);
    
    con = con(:);
    aon = aon(:);
    
end