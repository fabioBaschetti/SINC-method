function [conP,aonP,put,call] = SINC_DFTall(S,t,K,IR,DY,model,params,xi,j,Xc,Xm,N)

    if isrow(K)
        % do nothing
    else
        K = K';
    end
    
    Xl  = -Xc+Xm;
    Xh  =  Xc+Xm;
    
    wid =  Xh-Xl;           
    
    n   = 1:2:(N/2);            
    wn  = n/wid;            

    k   = log(K/S) - (IR-DY)*t - Xm;
    
    sn  = zeros(size(k,2),size(wn,2));
    cs  = zeros(size(k,2),size(wn,2));
    
    for i=1:size(k,2)
        sn(i,:) = sin(2*pi*k(i)*wn);
        cs(i,:) = cos(2*pi*k(i)*wn);
    end
    
    f1  = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn   ,t); 
    f2  = exp(-1i*2*pi*wn*Xm) .* charfun(model,params,xi,j,2*pi*wn-1i,t);
    
    ad1 = zeros(1,size(k,2));
    ad2 = zeros(1,size(k,2));
    
    for i=1:size(k,2)
        ad1(i) = sum((sn(i,:).*real(f1) - cs(i,:).*imag(f1)) ./ n);
        ad2(i) = sum((sn(i,:).*real(f2) - cs(i,:).*imag(f2)) ./ n);
    end
    
    conP = K*exp(-IR*t).*(0.5+(2/pi)*ad1);
    aonP = S*exp(-DY*t).*(0.5+(2/pi)*ad2);
    
    put = conP - aonP;  
    
    call = put+S*exp(-DY*t)-K*exp(-IR*t);
    
    conP = conP(:);
    aonP = aonP(:);
    put  = put(:);
    call = call(:);
    
end