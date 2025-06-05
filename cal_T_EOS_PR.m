function T = cal_T_EOS_PR(p,rho)

    maxStep = 100;
    tol = 1e-4;
    T0 = 300;
for i = 1:maxStep
    F = p - EOS_PR(rho,T0);
    f =  -EOS_PR_dT(rho,T0);
    delta = F ./ f;
    T_new = T0 - 0.9*delta;
    
    T0 = T_new; 
   if max(abs(delta)) < tol
        break;
    end
end

if i == maxStep
    warning('[cal_T_EOS_PR] Convergence tolerance not met. Results may be inaccurate.');
end   

    T = T0;
end