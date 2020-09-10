function [con,aon,put] = SINC_discFT_sbs(S,t,K,params_rHeston,j,xi,Xc,Xm,N)

    % OUTPUT: CoN, AoN and Plain Vanilla put prices from the FT-version
    % of the SINC formula. This uses as many evaluations of the CF as  
    % they are needed to reach a specified accuracy on CoN and AoN prices
    % (or, at most N if this is not met)

    if iscolumn(K)
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

    % initializing all quantities
    sn(:,1) = sin(2*pi*k*wn(1));
    cs(:,1) = cos(2*pi*k*wn(1));

    f1(1) = exp(-1i*2*pi*wn(1)*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*wn(1),t);                 
    f2(1) = exp(-1i*2*pi*wn(1)*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*(wn(1)-1i/(2*pi)),t);

    for m = 1:length(k)
        ad1(m,1) = sum((sn(m,1).*real(f1) - cs(m,1).*imag(f1)) ./ n(1));
        ad2(m,1) = sum((sn(m,1).*real(f2) - cs(m,1).*imag(f2)) ./ n(1));
    end

    con(:,1) = (0.5+2/pi*ad1(:,1));
    aon(:,1) = (0.5+2/pi*ad2(:,1));

    put(:,1) = K.*con(:,1) - S*aon(:,1);

    % tolerance on CoN and AoN prices
    tol = 1e-13;

    for i = 2:length(wn)
        
        % compute sin and cos in formula (10)
        sn(:,i) = sin(2*pi*k*wn(i));
        cs(:,i) = cos(2*pi*k*wn(i));
        
        % CF at points wn  for the CoN 
        f1(i) = exp(-1i*2*pi*wn(i)*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*wn(i),t);   
        % and wn-1i/(2*pi) for the AoN (exponential shift included)
        f2(i) = exp(-1i*2*pi*wn(i)*Xm) .* phirHeston(params_rHeston,xi,j,2*pi*(wn(i)-1i/(2*pi)),t);

        % compute sums in eq. (8) and (9) for CoN and Aon, resp.ly 
        for m = 1:length(k)
            ad1(m,i) = sum((sn(m,:).*real(f1) - cs(m,:).*imag(f1)) ./ n(1:i));
            ad2(m,i) = sum((sn(m,:).*real(f2) - cs(m,:).*imag(f2)) ./ n(1:i));
        end

        % conclude CoN and AoN put prices (eq. (8) and (9))
        con(:,i) = (0.5+2/pi*ad1(:,i));
        aon(:,i) = (0.5+2/pi*ad2(:,i));

        % put price as a linear combination of digital options
        put(:,i) = K.*con(:,i) - S*aon(:,i);

        % stopping condition (tolerance on CoN and AoN is met)
        if (mean(abs(con(:,i)-con(:,i-1)))<tol) && (mean(abs(aon(:,i)-aon(:,i-1)))<tol) 
            break
        end

    end

    con = con(:,end);
    aon = aon(:,end);
    put = put(:,end);

end