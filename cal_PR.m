%Calculate the coefficients of the PR equation and the temperature-related derivatives
%The medium is carbon dioxide (CO2)

function [a,b,R,dadT,d2adT2,cof] = cal_PR(T)
%%
Ru = 8.31443;
% CO2
MW = 44.0e-3;
Tc = 304.1282;
pc  = 7.3773e6;
omega  = 0.228;

k = 0.37464 + 1.54226*omega - 0.26992*omega^2;

R = Ru/MW;
a = 0.457235528921382*(R*Tc)^2/pc*(1+k*(1-sqrt(T/Tc))).^2;
b = 0.077796073903888*R*Tc./pc;
G = k*sqrt(T/Tc)./(1+k*(1-sqrt(T/Tc)));
dadT = -1./T.*a.*G;
d2adT2 =  0.457235528921382*R^2./T/2*k*(1+k)*Tc/pc.*sqrt(Tc./T);

%NASA 7 coefficient polynomials for CO2
cof = [2.35677352E+00, 8.98459677E-03,-7.12356269E-06, 2.45919022E-09 ,...
     -1.43699548E-13,-4.83719697E+04,9.90105222E+00]; 
end

