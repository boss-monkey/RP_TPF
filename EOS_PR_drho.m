function dpdrho = EOS_PR_drho(rho,T)

[a,b,R] = cal_PR(T);

dpdrho = R.*T./(1-b*rho).^2 - 2*a.*rho.*(1+b.*rho)./(1+2*b*rho-b^2*rho.^2).^2; 





