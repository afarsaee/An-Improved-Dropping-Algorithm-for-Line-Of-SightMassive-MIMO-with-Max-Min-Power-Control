function [H_dropped, n_user_dropped] = Drop_user_MRT(H, Ptot, channel_norm_input, n_max_drop, threshold_precision_sum_rate, flag_Ashkan1_Marzetta0)
[n_user,M_ant] = size(H);
HHH = find_rho_ij(H);
n_user_new = n_user;
flag_correlated_user = 1;
channel_norm_input_ref = channel_norm_input;
   while flag_correlated_user == 1 && n_user_new >= 2 && n_user_new > n_user - n_max_drop
        % find the correlated users indexes
        [cor_val,index_cor] = max(HHH(:));
        [cor_user1,cor_user2] = ind2sub(size(HHH),index_cor);
        HHH_new = (abs(HHH));
        HHH_new(cor_user1,cor_user2) = 0;
        HHH_new(cor_user2,cor_user1) = 0;
        channel_norm_sorted = zeros(1,n_user_new);
        index_write_norm = 1;
        for i = 1:n_user_new
            if i ~= cor_user1 && i ~= cor_user2
                channel_norm_sorted(index_write_norm) = channel_norm_input_ref(i);
                index_write_norm = index_write_norm + 1;
            end
        end
        channel_norm_sorted(n_user_new-1) = max([channel_norm_input_ref(cor_user1),channel_norm_input_ref(cor_user2)]);
        channel_norm_sorted(n_user_new)   = min([channel_norm_input_ref(cor_user1),channel_norm_input_ref(cor_user2)]);

        rho_threshold_ZF_current_SNR = find_threshold_CB_MaxMin_ComLetter_General(Ptot,n_user_new,channel_norm_sorted,threshold_precision_sum_rate);
        if cor_val > rho_threshold_ZF_current_SNR
            flag_drop = 1;
        else
            flag_drop = 0;
        end
        
        if flag_drop == 1
            diag_HHH = diag(abs(H*H'));
            HHH_user_1 = 1/sqrt(diag_HHH(cor_user1));
            HHH_user_2 = 1/sqrt(diag_HHH(cor_user2));
            if flag_Ashkan1_Marzetta0 == 1
                % Similar to Our paper
                compare_cor_user_1 = HHH_user_1*sum(HHH_new(cor_user1,:));
                compare_cor_user_2 = HHH_user_2*sum(HHH_new(cor_user2,:));
            elseif flag_Ashkan1_Marzetta0 == 0
                compare_cor_user_1 = max(HHH_new(cor_user1,:));
                compare_cor_user_2 = max(HHH_new(cor_user2,:));
            end
            % the new channel matrix
            H_new = zeros(n_user_new-1,M_ant);
            
            % keep the user which is less correlated to the others
            if compare_cor_user_1 > compare_cor_user_2
                H_new(1,:)      = H(cor_user2,:);
            else
                H_new(1,:)      = H(cor_user1,:);
            end
            
            % rebuild the channel matrix
            index_new = 1;
            for i = 2:n_user_new-1
                while index_new == cor_user1 || index_new == cor_user2
                    index_new = index_new + 1;
                end
                H_new(i,:) = H(index_new,:);
                index_new = index_new + 1;
            end
            
            % update the number of user and channel
            n_user_new = n_user_new - 1;
            H = H_new;
            HHH = find_rho_ij(H);
            channel_norm_input_ref = sqrt(diag(H*H'))';
       else
            flag_correlated_user = 0;
       end
   end
   n_user_dropped = n_user_new;   
   H_dropped      = H;
end