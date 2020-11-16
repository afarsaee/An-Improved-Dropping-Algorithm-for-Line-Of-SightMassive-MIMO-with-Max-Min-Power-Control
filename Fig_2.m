clc
clear all
close all
%%
addpath('func');
flag_write = 1;
gig = 10^9; c = 0.3 * gig; f0 = 1.9*gig; lambda = c/f0; d = lambda/2;
M = 64; n_user_ref = 3;
mySNR = 3;
alpha_LOS = 1;
theta_min = 30; theta_max = 150;
min_spacing_distance_user = lambda;
flag_correlation = 1;
R_min = 2000; R_max = 2000;
phi = exp(-4:0.01:3);
n_phi = length(phi);
sum_rate_ZF_not_dropped_maxmin = zeros(1,n_phi);
sum_rate_CB_not_dropped_maxmin = zeros(1,n_phi);
rho = zeros(1,n_phi);
Ptot = mySNR * n_user_ref;
%%
n_user_bound = 7;
Pn = 1;
alpha_bound = (n_user_bound/(n_user_bound+1))* log2(1+Ptot/(n_user_bound*Pn));
beta = 2/((Ptot/((2^alpha_bound)-1))-(n_user_bound-1));
%%
spacing_array = 1/2;
shadowing = 0;
for i_phi = 1:n_phi
    n_user = n_user_ref;
    separation_angle = phi(i_phi);
[channel_unit_norm,channel_gain,path_loss_dB] = gen_LOS_Channel(M,n_user,f0,R_max,R_min,...
                    alpha_LOS,theta_min, theta_max, min_spacing_distance_user, spacing_array,shadowing,flag_correlation,separation_angle);  
   H = channel_unit_norm';
   HHH = abs(H*H') - eye(n_user);
   rho(i_phi) = abs(HHH(1,2));
   %% ZF
   Hnormal = channel_unit_norm';
   UZFnormal_non_normalized = Hnormal'/(Hnormal*Hnormal');
   UZFnormal = zeros(size(UZFnormal_non_normalized));
   for i = 1:n_user_ref
      UZFnormal(:,i) = UZFnormal_non_normalized(:,i)/norm(UZFnormal_non_normalized(:,i));
   end
   GammaNormal = abs(diag(Hnormal*UZFnormal));
   sum_one_overGamma2Normal = sum(1./(GammaNormal.^2));
   for i = 1:n_user_ref
      sum_rate_ZF_not_dropped_maxmin(i_phi) = sum_rate_ZF_not_dropped_maxmin(i_phi) + log2(1+(Ptot/(sum_one_overGamma2Normal)));      
   end
   SINR_k_maxmin = myCB_MAXMIN(n_user_ref,Hnormal,Ptot,1e-5,1e-2);
   sum_rate_CB_not_dropped_maxmin(i_phi) = sum(log2(1+SINR_k_maxmin));
end
%%
[~,index] = sort(sum_rate_ZF_not_dropped_maxmin);
figure;
plot(rho(index),sum_rate_ZF_not_dropped_maxmin(index));
hold on
plot(rho(index),sum_rate_CB_not_dropped_maxmin(index));
legend('ZF','CB');
%%
if flag_write == 1
    for i_dummy = 1:1
            fZFnormal  = fopen('ZF_normal.txt','w');
            fCBnormal  = fopen('CB_normal.txt','w');
            n_write = length(sum_rate_ZF_not_dropped_maxmin);
            for i = 1:n_write
               fprintf(fZFnormal,'%0.6f %2.6f\n',rho(index(i))  ,real(sum_rate_ZF_not_dropped_maxmin(index(i))));
               fprintf(fCBnormal,'%0.6f %2.6f\n',rho(index(i))  ,real(sum_rate_CB_not_dropped_maxmin(index(i))));
            end
            fclose(fZFnormal);
            fclose(fCBnormal);
    end
end