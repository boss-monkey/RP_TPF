function p = EOS_PR(rho,T)

[a,b,R] = cal_PR(T);

p = rho.*R.*T./(1-b.*rho) - rho.^2.*a./(1+2*rho.*b-(b.*rho).^2);






