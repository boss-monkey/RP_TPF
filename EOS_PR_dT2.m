function d2p = EOS_PR_dT2(rho,T)

[a,b,R,dadT,d2adT2] = cal_PR(T);

d2p = -rho.^2.*d2adT2./(1+2*rho.*b-(b.*rho).^2);
end