function out = BS_closedForm(S,K,t,IR,DY,sigma,CP,type)

    if isrow(K)
        % do nothing
    else
        K = K';
        end

    x   = log(S./K) + (IR-DY).*t;
    sig = sigma .* sqrt(t);
    d1  = x./sig + sig./2;
    d2  = d1 - sig;
    pv  = exp(-IR.*t);

    put  = -S*exp(-DY*t).*normcdf(-d1) + pv.*K.*normcdf(-d2);
    call =  S*exp(-DY*t).*normcdf( d1) - pv.*K.*normcdf( d2);
    
    conP = pv.*K.*normcdf(-d2);
    conC = pv.*K.*normcdf( d2);
    
    aonP = exp(-DY*t)*S*normcdf(-d1);
    aonC = exp(-DY*t)*S*normcdf( d1);
    
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
