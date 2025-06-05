% Calculate energy
% Use NASA 7 coefficient polynomials as reference values for CO2
function e = cal_energy(rho,T)

[a,b,R,dadT,d2adT2,cof] = cal_PR(T);

    h_ideal = R.*T.*(cof(1) + cof(2)*T/2 + cof(3)*T.^2/3 + cof(4)*T.^3/4 + cof(5)*T.^4/5 + cof(6)./T);
    v = 1./rho;
    K1 = 1/sqrt(8)./b.*log((v+(1-sqrt(2)).*b)./(v+(1+sqrt(2)).*b));
    e = h_ideal - R.*T + (a - T.*dadT).*K1;


    end