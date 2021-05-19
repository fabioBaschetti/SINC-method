function x = frFFT_circ(p,e)
 
    if (e > 1) || (e < 0)
        error('eps does not belong in (0,1)')
    end
    
    N = length(p);
    
    eta = e / N;
    
    M = 0;
    while 2^M < 2*N
        M = M+1;
    end
    Q = 2^M - 2*N;
    
    q = zeros(1,2^M);
    q(1:N) = p .* exp(-1i*pi*eta*(-N/2:N/2-1).^2);
    
    C = [0:N-1 zeros(1,2^M-N)];
    for n = N+Q+2:2^M
        C(n) = -(2*N+Q)+n-1;
    end

    c = exp(1i*pi*eta*(C.^2));
    
    q_hat = fft(q);
    c_hat = fft(c);
    
    f_hat = ifft(q_hat .* c_hat);
    
    x = f_hat(1:N) .* exp(-1i*pi*eta*(-N/2:N/2-1).^2);
    
end
