function [rm,rp,p1,p2,p3,q1,q2,q3] = dh_Pade33_coeff(alpha,rho,a)
    
    aa  = sqrt(a.*(a+(0+1i)) - rho^2.*a.^2);
    rm  = -(0+1i)*rho.*a - aa;
    rp  = -(0+1i)*rho.*a + aa;
    
    b1  = -a .* (a+1i) / (2*gamma(1+alpha));
    b2  = (1-a*1i) .* a.^2 * rho / (2*gamma(1+2*alpha));  
    b3  = gamma(1+2*alpha) / gamma(1+3*alpha) * (a.^2.*(1i+a).^2 / (8*gamma(1+alpha)^2) + (a+1i).*a.^3*rho^2 / (2*gamma(1+2*alpha)));
    
    g0  = rm;
    g1  = -rm ./ (aa*gamma(1-alpha));
    g2  = rm ./ aa.^2 / gamma(1-2*alpha) .* (1 + rm./(2*aa)*gamma(1-2*alpha)/gamma(1-alpha)^2);
  
    den = g0.^3 + 2*b1.*g0.*g1 - b2.*g1.^2 + b1.^2.*g2 + b2.*g0.*g2;
    
    p1  = b1;
    p2  = (b1.^2.*g0.^2 + b2.*g0.^3 + b1.^3.*g1 + b1.*b2.*g0.*g1 - b2.^2.*g1.^2 + b1.*b3.*g1.^2 + b2.^2.*g0.*g2 - b1.*b3.*g0.*g2) ./ den;
    q1  = (b1.*g0.^2 + b1.^2.*g1 - b2.*g0.*g1 + b3.*g1.^2 - b1.*b2.*g2 - b3.*g0.*g2) ./ den;
    q2  = (b1.^2.*g0 + b2.*g0.^2 - b1.*b2.*g1 - b3.*g0.*g1 + b2.^2.*g2 - b1.*b3.*g2) ./ den;
    q3  = (b1.^3 + 2*b1.*b2.*g0 + b3.*g0.^2 - b2.^2.*g1 + b1.*b3.*g1) ./ den;
    p3  = g0 .* q3;

end