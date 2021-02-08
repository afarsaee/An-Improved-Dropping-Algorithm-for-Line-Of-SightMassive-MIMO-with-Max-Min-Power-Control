function thr_CB_maxmin = find_threshold_CB_MaxMin_ComLetter_General(Ptot,K,channel_norm,threshold_precision_sum_rate)
    %% This function finds the threshold value for CB
% The threshold specifies the region in which it is better
% to drop the correlated user from service
% if measured_rho > "rho" --> R_dropped > R_no_drop
% if measured_rho < "rho" --> R_dropped < R_no_drop

% The equalized SINR for each user for a correlated channel 
% of size K is (user K and K-1 are correlated):
% SINR_i=\frac{\beta \sigma_K^2}{\beta \rho^2 \sigma_K^2 + N0}

% In this code, we use bisection methoed to search for the
% optimal threshold, given K,Ptot,N0,beta,
% so we find R_no_drop for each guess of "rho"
% then try to find where "|R_no_drop - R_drop| < threshold_precision"
% which is to avoid a large number of steps (less accurate)
N0 = 1;
flag_run = 1;
rho2 = 0.5; % first guess

sum_one_over_norm_user_1_to_K_minus_2 = sum( (1./channel_norm(1:K-2)).^2 );
R_DU_denom_part1 = sum_one_over_norm_user_1_to_K_minus_2;
% Ideally, in FP, the user with the lowest channel norm is dropped:
R_DU_denom_part2 = 1/((channel_norm(K-1))^2);
R_dropped_ref = (K-1)*log2(1+Ptot/(R_DU_denom_part1+R_DU_denom_part2));

rho0 = 0;
rho1 = 1;
%% Main Loop
while flag_run == 1
        gamma1 = sum_one_over_norm_user_1_to_K_minus_2;
        gamma2  = (1/((channel_norm(K-1))^2)) + (1/((channel_norm(K))^2));
        b_z = -(Ptot * rho2 + N0 *(gamma1 + gamma2));
        c_z = Ptot;
        a_z = N0 * gamma1 * rho2;
        
        z = (-b_z - sqrt(b_z^2 - 4*a_z*c_z))/(2*a_z);
        R_iteration = K*log2(1 + z);
    
    if R_iteration < R_dropped_ref
        rho1 = rho2;
        rho2 = (rho0 + rho2)/2;
    else
        rho0 = rho2;
        rho2 = (rho1+rho2)/2;
    end
    if abs(R_iteration - R_dropped_ref) < threshold_precision_sum_rate || abs(rho0 - rho1) < 1e-2
        flag_run = 0;
        break;
    end
    if rho2 < threshold_precision_sum_rate
        rho2 = 0;
        flag_run = 0;
        break;
    end
    if rho2 > 1-threshold_precision_sum_rate
        rho2 = 1;
        flag_run = 0;
        break;
    end
end
thr_CB_maxmin = sqrt(rho2);
end