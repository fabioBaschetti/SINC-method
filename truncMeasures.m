function [Xc,Xm] = truncMeasures(params_rHeston,j,xi,t,bl,bh,N,L)

    % OUTPUT: values Xc and Xm to determine the truncation range [Xl,Xh] 
    % used with SINC formulas
 
    % initial value for Xc based on guesses bl and bh 
    % of Xl and Xh respectively
    bc = (bh-bl)/2;
        
    % sampled Fourier modes \omega_n for the moments
    n   = 1:N;
    wn  = n/(2*bc);
    
    % CF at points wn 
    phi = phirHeston(params_rHeston,xi,j,2*pi*wn,t);
    
    % numerical moments as of Appendix E
    m1  = -2*bc * sum(imag(phi) .* (-1).^n ./ (n*pi));
    m2  = bc^2/3 + 4*bc^2 * sum(real(phi) .* (-1).^n ./ (n*pi).^2);
    m3  = -2*bc^3 * sum(imag(phi) .* (-1).^n ./ (n*pi) .* (1-6./(n*pi).^2));
    m4  = bc^4/5 + 8*bc^4 * sum(real(phi) .* (-1).^n ./ (n*pi).^2 .* (1-6./(n*pi).^2));
    
    % from moments to cumulants
    c1  = m1;
    c2  = m2 - m1^2;
    c4  = m4 - 4*m1*m3 + 6*m1^2*m2 - 3*m1^4;
    
    % Fang-Oosterlee truncation rule
    Xl  = c1 - L*sqrt(c2+sqrt(c4));
    Xh  = c1 + L*sqrt(c2+sqrt(c4));
    
    % compute Xc and Xm to be used in the SINC
    Xc  = (Xh - Xl)/2;
    Xm  = (Xh + Xl)/2;

end