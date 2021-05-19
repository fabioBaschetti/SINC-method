function [x,pdf] = chf2pdf(model,params,t,xi,j,Xc,Xm,N)

    cf = @(u) charfun(model,params,xi,j,u,t);
    
    Xl  = -Xc + Xm;
    Xh  =  Xc + Xm;

    wid =  Xh - Xl;
    
    n   = -N/2:N/2-1;
    wn  = n / wid;
    
    phi = cf(2*pi*wn).*exp(-1i*2*pi*wn*Xm);
    phi(N/2+1) = 1;
    phi = phi/wid;

    pdf = fft(phi);
    pdf = pdf(1:2:end);
    M   = size(pdf,2);
    pdf = [pdf(M/2+1:M) pdf(1:M/2)];
    pdf = real(pdf);

    dx  = wid/M;
    x   = (-M/2:(M/2)-1)*dx + Xm;
        
end
