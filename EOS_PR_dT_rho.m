function d2p = EOS_PR_dT_rho(rho,T)

[a,b,R,dadT,d2adT2] = cal_PR(T);

d2p = R./(1-b*rho).^2 - 2*dadT.*rho.*(1+b.*rho)./(1+2*b*rho-b^2*rho.^2).^2; 

end