function out = BSImpliedVolPut(S,K,T,r,P)

    n      =  size(K,2);
    
    sigmaL = 1e-10;
    PL     = BSFormulaPut(S,K,T,r,sigmaL);

    sigmaH = 10;
    PH     = BSFormulaPut(S,K,T,r,sigmaH);
    
    while mean(sigmaH - sigmaL) > 1e-10
        
        sigma  = (sigmaL+sigmaH) ./ 2;
        PM     = BSFormulaPut(S,K,T,r,sigma);
        
        PL     = PL + (PM < P) .* (PM-PL);
        sigmaL = sigmaL + (PM < P) .* (sigma-sigmaL);
        
        PH     = PH + (PM >= P).* (PM-PH);
        sigmaH = sigmaH + (PM >= P).* (sigma-sigmaH);
        
    end
    
    out = sigma;
    
end

function out = BSFormulaPut(S,K,T,r,sigma)

     x   = log(S./K) + r.*T;
     sig = sigma .* sqrt(T);
     d1  = x./sig + sig./2;
     d2  = d1 - sig;
     pv  = exp(-r.*T);
     
     put = pv.*K.*normcdf(-d2) - S.*normcdf(-d1);

     out = put;

end