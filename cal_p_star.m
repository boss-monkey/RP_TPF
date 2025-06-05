%This program is figure 4 in the paper.

function star = cal_p_star(W_R,W_L,maxStep,tol)
if nargin < 3
    maxStep = 100;
    tol = 10;
end

p_star = (W_L(3,:).*W_R(3,:)).^0.5; %Initial value
rho0 = (W_L(1,:).*W_R(1,:)).^0.5;   %initial value

rho_L_split = cal_split(rho0, W_L(1,:)); %split density initials to enhance computational stability
rho0_L_plus = rho_L_split(1,:);
rho0_L_minus =  rho_L_split(2,:);

rho_R_split = cal_split(rho0, W_R(1,:));
rho0_R_plus = rho_R_split(1,:);
rho0_R_minus =  rho_R_split(2,:);

for i = 1:maxStep
    u_R   = W_R(2,:);
    u_L   = W_L(2,:);

    p_L_split = cal_split(p_star, W_L(3,:));  %pressure splitting strategy
    p_star_L_plus = p_L_split(1,:);
    p_star_L_minus =  p_L_split(2,:);
    
    p_R_split = cal_split(p_star, W_R(3,:));
    p_star_R_plus = p_R_split(1,:);
    p_star_R_minus =  p_R_split(2,:);
    
    
    F_L = cal_F(W_L, p_star_L_plus,rho0_L_plus);
    J_L =  cal_J(W_L, p_star_L_minus,rho0_L_minus);

    F_R = cal_F(W_R, p_star_R_plus,rho0_R_plus);
    J_R = cal_J(W_R, p_star_R_minus,rho0_R_minus);

    dF_L = cal_dFdp(W_L,p_star_L_plus,rho0_L_plus);
    dJ_L = cal_dJdp(W_L, p_star_L_minus,rho0_L_minus);

    dF_R = cal_dFdp(W_R, p_star_R_plus,rho0_R_plus);
    dJ_R = cal_dJdp(W_R, p_star_R_minus,rho0_R_minus);

%%    
    G = (u_L-u_R)+(F_L+J_L)+(F_R+J_R);
    dGdp = dF_L+dJ_L+dF_R+dJ_R;
    delta = G./(dGdp);
    p_star_new = p_star - delta;
    
    p_star  = p_star_new;
   
   if max(abs(delta)) < tol
        break;
   end
  
end
disp(['Calculate the intermediate pressure of the RP, current absolute error ' num2str(max(abs(delta)))  '，max error' num2str(tol)]); 
star = real([p_star;F_L;F_R;J_L;J_R]);

end
