function dpdT = EOS_PR_dT(rho,T)

[a,b,R,dadT] = cal_PR(T);

dpdT = rho.*R./(1-b.*rho) - rho.^2.*dadT./(1+2*rho.*b-(b.*rho).^2);







