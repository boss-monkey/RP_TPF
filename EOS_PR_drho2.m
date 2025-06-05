function d2p = EOS_PR_drho2(rho,T)

[a,b,R] = cal_PR(T);

d2p = 2.*b.*R.*T./(1-b*rho).^3 - 2*a.*(1+3*b.^2 .*rho.^2+2*b.^3.*rho.^3)...
    ./(1+2*b*rho-b^2*rho.^2).^3;
end