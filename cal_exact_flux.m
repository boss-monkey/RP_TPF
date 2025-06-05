%This program calculates the flow variables on each interface.
function flux = cal_exact_flux(W_L,W_R,star,x_ref)

p_star = star(1,:);
F_L = star(2,:); F_R = star(3,:);
J_L = star(4,:); J_R = star(5,:);

rho_L = W_L(1,:); rho_R = W_R(1,:);
u_L = W_L(2,:); u_R = W_R(2,:);
p_L= W_L(3,:); p_R= W_R(3,:);
T_L = W_L(4,:); T_R = W_R(4,:);

rho0 = (W_L(1,:).*W_R(1,:)).^0.5;


    if abs(p_star-p_L)+abs(p_star-p_R) <=1e-9
        condition = 1;   %No shock waves or rarefaction waves
    elseif (p_star-p_L)>0 && (p_star-p_R)> 0
        condition = 2;   %Double shock wave
    elseif (p_star-p_L)<0 && (p_star-p_R)< 0
        condition = 3;   %Double rarefaction wave
    elseif (p_star-p_L)<0 && (p_star-p_R)> 0
        condition = 4;    %Left rarefaction wave, right shock wave
    elseif (p_star-p_L)>0 && (p_star-p_R)< 0
        condition = 5;    %Left shock wave, right expansion wave
    end

        if condition == 1
            S_L_head = 0;  
            S_L_tail = 0;
            S_R_head = 0;
            S_R_tail = 0;

        elseif condition == 2

            rho_star_L = cal_rho_hugoniot(W_L, rho0, p_star);
            rho_star_R = cal_rho_hugoniot(W_R, rho0, p_star);
            S_L_head = cal_shock_speed(W_L, rho0,1, p_star);
            S_R_head = cal_shock_speed(W_R, rho0,-1, p_star);
            u_star = ((u_L+F_L)+(u_R-F_R))/2;       %Take the average to reduce calculation errors.
            T_star_L = cal_T_EOS_PR(p_star,rho_star_L); 
            T_star_R = cal_T_EOS_PR(p_star,rho_star_R);


            S_L_tail = S_L_head;
            S_R_tail = S_R_head;
 
        elseif condition == 3
            s_L = cal_entropy(rho_L,T_L);
            s_R = cal_entropy(rho_R,T_R);
            rho_star_L =  cal_rho_isentropic(s_L, rho0, p_star);
            rho_star_R = cal_rho_isentropic(s_R, rho0, p_star);
            u_star = ((u_L+J_L)+(u_R-J_R))/2;       
            T_star_L = cal_T_EOS_PR(p_star,rho_star_L);
            T_star_R = cal_T_EOS_PR(p_star,rho_star_R);
            c_star_L = cal_c_rho_T(rho_star_L,T_star_L);
            c_star_R = cal_c_rho_T(rho_star_R,T_star_R);
            c_L = cal_c_rho_T(rho_L,T_L);
            c_R = cal_c_rho_T(rho_R,T_R);

            S_L_head = u_L - c_L;
            S_L_tail = u_star - c_star_L;
            S_R_head = u_R + c_R;
            S_R_tail = u_star + c_star_R;

        elseif condition == 4   
            s_L = cal_entropy(rho_L,T_L); 
            rho_star_L = cal_rho_isentropic(s_L, rho0, p_star);
            rho_star_R = cal_rho_hugoniot(W_R, rho0, p_star);
            u_star = ((u_L+J_L)+(u_R-F_R))/2;      
            T_star_L = cal_T_EOS_PR(p_star,rho_star_L); 
            T_star_R = cal_T_EOS_PR(p_star,rho_star_R);

            c_star_L = cal_c_rho_T(rho_star_L,T_star_L);
            c_L = cal_c_rho_T(rho_L,T_L);

            S_L_head = u_L - c_L;
            S_L_tail = u_star - c_star_L;
            S_R_head = cal_shock_speed(W_R, rho0,-1, p_star);
            S_R_tail = S_R_head;

        elseif condition == 5
            rho_star_L = cal_rho_hugoniot(W_L, rho0, p_star);
            s_R = cal_entropy(rho_R, T_R);
            rho_star_R = cal_rho_isentropic(s_R, rho0, p_star);
            u_star = ((u_L+F_L)+(u_R-J_R))/2;
            T_star_L = cal_T_EOS_PR(p_star,rho_star_L); 
            T_star_R = cal_T_EOS_PR(p_star,rho_star_R);
            c_star_R = cal_c_rho_T(rho_star_R,T_star_R);
            c_R = cal_c_rho_T(rho_R,T_R);

            S_L_head = cal_shock_speed(W_L, rho0,1, p_star);
            S_L_tail = S_L_head;
            S_R_head = u_R + c_R;
            S_R_tail = u_star + c_star_R;

         end

        if S_L_head>=x_ref                          % Interface is located in the left undisturbed region
                 W = W_L;                           
        elseif S_L_head<x_ref && x_ref<S_L_tail     % Interface is located in the left rarefaction wave
                W = cal_wave_p(W_L,p_star,rho_star_L,1, x_ref);
        elseif S_L_tail<x_ref &&  x_ref<S_R_tail   
            if u_star >= x_ref                      % Interface is at the left contact discontinuity
                W = [rho_star_L;u_star;p_star;T_star_L];
            else                                    % Interface is at the right contact discontinuity
                W = [rho_star_R;u_star;p_star;T_star_R];
            end
        elseif S_R_head>x_ref && x_ref>S_R_tail     %Interface is located in the right rarefaction wave
                W = cal_wave_p(W_R,p_star,rho_star_R,-1, x_ref);
        elseif x_ref>=S_R_tail 
                W = W_R;                           %Interface is located in the right undisturbed region
        end
        flux = W;

end