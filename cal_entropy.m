% Calculate entropy
% Use NASA 7 coefficient polynomials as reference values for CO2
function s = cal_entropy(rho,T)

[a,b,R,dadT,d2adT2, cof] = cal_PR(T);

    s_ideal = R.*(cof(1)*log(T) + cof(2)*T + cof(3)*T.^2/2 + cof(4)*T.^3/3 + cof(5)*T.^4/4 + cof(7));
    v = 1./rho;
    K1 = 1/sqrt(8)./b.*log((v+(1-sqrt(2)).*b)./(v+(1+sqrt(2)).*b));
    s = s_ideal + R.*log((v-b)./b) - K1.*dadT;

    end