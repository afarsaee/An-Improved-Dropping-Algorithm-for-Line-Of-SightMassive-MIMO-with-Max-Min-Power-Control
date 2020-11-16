clc
clear all
close all
%% This code can be used to run a simulation as in Fig. 3 and 4 of our paper
% This code's scenario is slightly different from our paper
%% Flags in code
flag_write = 1;
%% Initialization
n_channel = 10000;
rng('default');
addpath('func');                % adding the path for func
%% LOS Config.
M = 40;
n_user_ref = 8;
gig = 10^9; c = 0.3 * gig; f0 = 30*gig; lambda = c/f0;
phi_min = 30; phi_max = 150;
R_min = 200; R_max = 200;% Simulation parameters,
min_spacing_x = 0.01;
mySNRdB = 0:5:25;
index_show = find(mySNRdB == 25);
mySNR = 10.^(mySNRdB/10);
n_SNR = length(mySNR);
thr_Marzetta = 0.75*ones(1,n_SNR);
%% Variables
% sum-rate variables
sum_rate_ZF              = zeros(1,n_SNR);
sum_rate_ZF_Proposed     = zeros(1,n_SNR);
CDFsum_rate_ZF           = zeros(n_channel,n_SNR);
CDFsum_rate_ZF_Proposed  = zeros(n_channel,n_SNR);
% using lambda spacing
spacing_array = 0.5;
% repeat the simulation for a same channel for all the SNRs
% so we need to record the "Href" once, then repeat the sims
% for it
Href                = zeros(n_channel,M,n_user_ref);
%% Main Loop
accum_error = zeros(n_channel,n_SNR);
wrong_drop = 0;
% record the number of dropped users
n_drop = zeros(n_channel,n_SNR);
for i_SNR = 1:n_SNR
    % finding corresponding total power
    Ptot = mySNR(i_SNR) * n_user_ref;
%% Repeat the simulation for n_channel realizations
    for i_channel = 1:n_channel
        % generate the channel for the first time
       if i_SNR == 1
            n_user = n_user_ref;
            [channel_unit_norm,~,path_loss_dB] = Gen_LOS_Channel_General(M,n_user,f0,R_max,R_min,...
                            phi_min, phi_max, min_spacing_x, spacing_array);
            Href(i_channel,:,:) = channel_unit_norm; % we read the normalized channel
            % if channel_unit_norm is used for the channel matrix
            % then, it is implicitly assumed that TX power at the BS is
            % scaled down when M increases to have a fair comparison to
            % small M values
       else
           % or read the same channel for the next SNR
           % and also read the channel gain
            channel_unit_norm = squeeze(Href(i_channel,:,:));
       end
       %% Read the channel matrix
       HDL_unit_norm      = channel_unit_norm';
       H = HDL_unit_norm;
       %% ZF Proposed
       [H_dropped,n_user_dropped] = Drop_user_ZF(H, Ptot, ones(1,n_user), 1, 5, 1);
       UZF_non_normalized = pinv(H_dropped);
       sum_filter_norm2 = sum(diag(UZF_non_normalized'*UZF_non_normalized));
       SNR_ZF_CDA = Ptot/(sum_filter_norm2);
       sum_rate_ZF_Proposed(i_SNR) = sum_rate_ZF_Proposed(i_SNR) + (n_user_dropped) * log2(1+SNR_ZF_CDA);
       CDFsum_rate_ZF_Proposed(i_channel,i_SNR) = (n_user_dropped) * log2(1+SNR_ZF_CDA);
       %% Normal one:
       n_user = n_user_ref;
       H = HDL_unit_norm;
       UZF_non_normalized = pinv(H);
       sum_filter_norm2 = sum(diag(UZF_non_normalized'*UZF_non_normalized));
       SNR_ZF = Ptot/(sum_filter_norm2);
       sum_rate_ZF(i_SNR) = sum_rate_ZF(i_SNR) + n_user_ref * log2(1+SNR_ZF);
       CDFsum_rate_ZF(i_channel,i_SNR) = n_user_ref * log2(1+SNR_ZF);
       %%
       if CDFsum_rate_ZF(i_channel,i_SNR) > CDFsum_rate_ZF_Proposed(i_channel,i_SNR)
           wrong_drop = wrong_drop + 1;
           accum_error(i_channel,i_SNR) = accum_error(i_channel,i_SNR) + (CDFsum_rate_ZF(i_channel,i_SNR)-CDFsum_rate_ZF_Proposed(i_channel,i_SNR));
       end
    end
end
sum_accum_error = sum(accum_error);
sum_n_drop      = sum(n_drop);
avg_accum_error = sum_accum_error./sum_n_drop;
sum_rate_ZF          = sum_rate_ZF/n_channel;
sum_rate_ZF_Proposed = sum_rate_ZF_Proposed/n_channel;   
%% plot the average sum-rate over SNR
figure;
plot(mySNRdB,sum_rate_ZF);
hold on;
plot(mySNRdB,sum_rate_ZF_Proposed);
legend('ZF normal','ZF proposed');
%% displaying the avg sum-rate
display(['avg proposed = ',num2str(sum_rate_ZF_Proposed(index_show))]);
%% Writing the results
if flag_write == 1
    for i_dummy = 1:1
            %%
            name_dropped     = sprintf('Fig_3_4_ZFsumrate_Proposed_%d_%d.txt',M,n_user_ref);
            name_not_dropped = sprintf('Fig_3_4_ZFsumrate_not_dropped_%d_%d.txt',M,n_user_ref);
            
            fsumratenotdropped          = fopen(name_not_dropped,'w');
            fsumratedropped             = fopen(name_dropped,'w');
            n_write = length(sum_rate_ZF_Proposed);
            for i = 1:n_write
               fprintf(fsumratenotdropped,'%0.6f %2.6f\n', mySNRdB(i) ,sum_rate_ZF(i));
               fprintf(fsumratedropped,'%0.6f %2.6f\n', mySNRdB(i) ,sum_rate_ZF_Proposed(i));
            end
            fclose(fsumratenotdropped);
            fclose(fsumratedropped);
    end
end