function out = dh_Pade33(y,rm,rp,p1,p2,p3,q1,q2,q3)

    h_pade = (p1*y + p2*y.^2 + p3*y.^3) ./ (1 + q1*y + q2*y.^2 + q3*y.^3);
    dh  = 1/2 * (h_pade-rm) .* (h_pade-rp);
    
    out = dh;

end