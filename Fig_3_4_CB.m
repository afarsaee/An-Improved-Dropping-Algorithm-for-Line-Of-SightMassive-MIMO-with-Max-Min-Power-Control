clc
clear all
close all
%% This code can be used to run a simulation as in Fig. 3 and 4 of our paper
% This code's scenario is slightly different from our paper
%% Flags
flag_write = 1;
%% Initialization
n_channel = 10000;
threshold_precision_sum_rate = 1e-4; threshold_precision_for_threshold = 1e-4;
Ptot_margin = 0.01; thr_error_SINR_CB = 1e-3;
rng('default');
addpath('func');                % adding the path for func
%% LOS Config.
M = 40; n_user_ref = 8;
gig = 10^9; c = 0.3 * gig; f0 = 30*gig; lambda = c/f0;
theta_min = 30; theta_max = 150;
spacing_array = 0.5; % using lambda spacing
min_spacing_x = 0.01;
R_min = 200; R_max = 200;
mySNRdB = 0:5:25;
mySNR = 10.^(mySNRdB/10);
n_SNR = length(mySNR);
%% Variables
sum_rate_CB             = zeros(1,n_SNR);
sum_rate_CB_Proposed     = zeros(1,n_SNR);
CDFsum_rate_CB          = zeros(n_channel,n_SNR);
CDFsum_rate_CB_Proposed  = zeros(n_channel,n_SNR);
Href = zeros(n_channel,M,n_user_ref);
%% Main Loop
N0 = 1;
wrong_drops = 0;
accum_error = zeros(n_channel,n_SNR);
for i_SNR = 1:n_SNR
    % finding corresponding total power
    % SNR_FP = beta*Ptot/n_user_ref*N0 --> (N0 = beta = 1)
    Ptot = mySNR(i_SNR) * n_user_ref;
%% Repeat the simulation for n_channel realizations
    for i_channel = 1:n_channel
        % generate the channel for the first time
       if i_SNR == 1
       n_user = n_user_ref;
       [channel_unit_norm,~,path_loss_dB] = Gen_LOS_Channel_General(M,n_user,f0,R_max,R_min,...
                            theta_min, theta_max, min_spacing_x, spacing_array);
            Href(i_channel,:,:) = channel_unit_norm;
       else
           % or read the same channel for the next SNR
            channel_unit_norm = squeeze(Href(i_channel,:,:));
       end
       %% Read the channel
       HDL_unit_norm = channel_unit_norm';
       H = HDL_unit_norm;
        
       HHH = find_rho_ij(H);
       [H_drop,n_user_drop] = Drop_user_MRT(H, Ptot, ones(1,n_user_ref), 5, threshold_precision_sum_rate, 1);
       %% CB Dropped
       [SINR_k_maxmin_CB_Proposed,~] = myCB_MAXMIN(n_user_drop,H_drop,Ptot,threshold_precision_sum_rate, Ptot_margin);
       if max(SINR_k_maxmin_CB_Proposed(:)) - min(SINR_k_maxmin_CB_Proposed(:)) > thr_error_SINR_CB
           error('error!');
       end
       sum_rate_CB_Proposed(i_SNR) = sum_rate_CB_Proposed(i_SNR) + sum(log2(1+SINR_k_maxmin_CB_Proposed));
       CDFsum_rate_CB_Proposed(i_channel,i_SNR) = sum(log2(1+SINR_k_maxmin_CB_Proposed));
       %% No Dropping
       n_user = n_user_ref;
       H = HDL_unit_norm;
       [SINR_k_maxmin_CB,~] = myCB_MAXMIN(n_user,H,Ptot,threshold_precision_sum_rate, Ptot_margin);
       % check if the SINRs are not equal!
       if max(SINR_k_maxmin_CB(:)) - min(SINR_k_maxmin_CB(:)) > thr_error_SINR_CB
           error('error!');
       end
       sum_rate_CB(i_SNR) = sum_rate_CB(i_SNR) + sum(log2(1+SINR_k_maxmin_CB));
       CDFsum_rate_CB(i_channel,i_SNR) = sum(log2(1+SINR_k_maxmin_CB));
       if n_user_drop == n_user_ref && abs(CDFsum_rate_CB(i_channel,i_SNR) - CDFsum_rate_CB_Proposed(i_channel,i_SNR)) > 0.1
           error('Error!');
       end
       if CDFsum_rate_CB_Proposed(i_channel,i_SNR) < CDFsum_rate_CB(i_channel,i_SNR)
           wrong_drops = wrong_drops + 1;
           accum_error(i_channel,i_SNR) = accum_error(i_channel,i_SNR) + (CDFsum_rate_CB(i_channel,i_SNR)-CDFsum_rate_CB_Proposed(i_channel,i_SNR));
       end
    end
end
accum_error_avg = mean(accum_error);
sum_rate_CB      = sum_rate_CB/n_channel;
sum_rate_CB_Proposed = sum_rate_CB_Proposed/n_channel;
%%
figure;
plot(mySNRdB,sum_rate_CB);
hold on;
plot(mySNRdB,sum_rate_CB_Proposed);
legend('CB normal','CB proposed');
%%
if flag_write == 1
    for i_dummy = 1:1
            %%
            name_dropped     = sprintf('Fig_5_CBsumrate_Proposed_%d_%d.txt',M,n_user_ref);           
            name_not_dropped = sprintf('Fig_5_CBsumrate_not_dropped_%d_%d.txt',M,n_user_ref);
            fsumratenotdropped          = fopen(name_not_dropped,'w');
            fsumratedropped             = fopen(name_dropped,'w');
            n_write = length(sum_rate_CB_Proposed);
            for i = 1:n_write
               fprintf(fsumratenotdropped,'%0.6f %2.6f\n', mySNRdB(i) ,sum_rate_CB(i));
               fprintf(fsumratedropped,'%0.6f %2.6f\n', mySNRdB(i) ,sum_rate_CB_Proposed(i));
            end
            fclose(fsumratenotdropped);
            fclose(fsumratedropped);
    end
end