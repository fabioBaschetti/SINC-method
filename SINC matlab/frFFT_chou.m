function out = frFFT_chou(h,e)

    if (e > 1) || (e < 0)
        error('eps does not belong in (0,1)')
    end

    N    = length(h);
    
    beta = e / N;
    
    y    = [exp(-1i*pi*(0:N-1).^2*beta).*h, zeros(1,N)];
    z    = [exp( 1i*pi*(0:N-1).^2*beta)   , exp( 1i*pi*(N:-1:1).^2*beta)];
    
    Dy   = fft(y);
    Dz   = fft(z);
    
    xx   = Dy.*Dz;
    
    ii   = ifft(xx);
    
    p    = [exp(-1i*pi*(0:N-1).^2*beta)   , zeros(1,N)];
    
    hh   = p.*ii;
    hh   = hh(1:N);
    
    out  = hh;
    
end