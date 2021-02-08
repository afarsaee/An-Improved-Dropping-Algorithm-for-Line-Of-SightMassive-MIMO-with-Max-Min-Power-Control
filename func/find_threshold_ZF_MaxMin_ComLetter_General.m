function thr_ZF_maxmin = find_threshold_ZF_MaxMin_ComLetter_General(Ptot,n_user,channel_norm,flag_cell_edge)
%% see (19) and (20) of the letter to check the threshold
% alpha = R_DU/K
% (1/||h_K||^2) + (1/||h_{K-1}||^2)
% \sum_{i=1}^{i=K-2} (1/||h_i||^2)
if flag_cell_edge == 1
    channel_norm = ones(1,n_user);
end
    sum_one_over_norm_user_1_to_K_minus_2 = sum( (1./channel_norm(1:n_user-2)).^2 );
    sum_one_over_norm_user_KKminus1 = sum((1./channel_norm(n_user-1:n_user)).^2);
    R_DU_denom_part1 = sum_one_over_norm_user_1_to_K_minus_2;
    R_DU_denom_part2 = 1/((channel_norm(n_user-1))^2);
    R_DU = (n_user-1)*log2(1+Ptot/(R_DU_denom_part1+R_DU_denom_part2));
    alpha_ZF_power = (R_DU/n_user);
    denom_ZF_threshold = (Ptot/(((2^alpha_ZF_power) - 1))) - sum_one_over_norm_user_1_to_K_minus_2;
    first_term = sum_one_over_norm_user_KKminus1/(denom_ZF_threshold);
    if 1-first_term > 0
        thr_ZF_maxmin = sqrt(1-first_term);
    else
        thr_ZF_maxmin = 0;
    end
end