clearvars
clc
%% Set initial values
p_0 = 1.5e+07;  p_1 = 1.5e+06;
rho_0 = 200;    rho_1 = 20;
u_0 =  50;        u_1 = 0;
T_0 = cal_T_EOS_PR(p_0,rho_0);
T_1 = cal_T_EOS_PR(p_1,rho_1);

%% Set calculation domain and calculation time
time = 5e-5;                    %time
x1 = 0; x2 = 0.1; x3 = 0.05;    %range
N = 400;                        %grid number
x = linspace(x1,x2,N);

%% Calculat the Intermediate pressure
W_L = [rho_0; u_0 ;p_0; T_0];
W_R = [rho_1; u_1 ;p_1; T_1];   %The left and right state
star = cal_p_star(W_R,W_L);

disp('Intermediate pressure calculation successful, start calculating interface flux'); 
%% Calculate the flow variables on each interface
W_1 = zeros(4,N);
x_ref =(x-x3)/time;
for i = 1:N
    flux = cal_exact_flux(W_L,W_R,star,x_ref(i));
    W_1(:,i) = flux;
    disp(['Current grid interface' num2str(i)]); 
end
%% Plot
rho1 = W_1(1,:); u1 = W_1(2,:);
p1   = W_1(3,:); T1   = W_1(4,:);
c1 = cal_c_rho_T(rho1,T1);
result = [x',rho1',u1',p1',T1'];

figure(1);
subplot(2,2,1)      % rho
plot(x,rho1,'-r','linewidth',1)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('x');
ylabel('rho');
subplot(2,2,2)      % u
plot(x,u1,'-b','linewidth',1)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('x');
ylabel('u');
subplot(2,2,3)      % pressure
plot(x,p1,'-m','linewidth',1)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('x');
ylabel('p');
subplot(2,2,4)     % temperature
plot(x,T1,'-k','linewidth',1)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('x');
ylabel('T');



