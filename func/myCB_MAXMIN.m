function [SINR_k_maxmin, iterations, SINR_k_equ, Ptot_consumed] = myCB_MAXMIN(n_user,h,Ptot,diff_compare_threshold, Ptot_margin)
% #iterations to find the max-min SINR
n_iteration = 100;

% normalized correlation matrix
hnormalized = zeros(size(h));
for i = 1:n_user
   hnormalized(i,:) = h(i,:)/norm(h(i,:));
end
rho2_ij = abs(hnormalized*hnormalized').^2;

% find the norm of the channel vectors
HHH = h*h';
norm_hi2 = abs(diag(HHH)); % ||h||_i^2
b1 = (1./norm_hi2); % 1./||h||^2

% the correlation matrix
RhoijMatrix = rho2_ij - diag(diag(rho2_ij));   

% the lower bound for the SINR
a_min = 0.01;                                   

% the upper bound for the SINR
a_max = (Ptot/n_user)*max(norm_hi2)*n_user;  

% identity matrix
IK  = eye(n_user);

SINR_candidate  = (a_min+a_max)/2;
max_min_power = SINR_candidate * ((IK-SINR_candidate*RhoijMatrix)\b1);
index = 0;
%% Main loop to find the max-min SINR and the max-min power control
% continue the loop till the #iterations and while the SINR accuracy is not met
while index <= n_iteration
   index = index + 1;
   
   % if all the power coefficients are positive and the sum of them is less
   % than the Ptot, then increase the interval of the search (we can
   % increase the power coefficients, right?)
   if sum(sign(max_min_power)) == n_user && sum(abs(max_min_power)) < Ptot
      a_min  = SINR_candidate;
      SINR_candidate = (a_min+a_max)/2;
   else
   % otherwise, decrease the search interval
      a_max  = SINR_candidate;
      SINR_candidate = (a_min+a_max)/2;
   end
   
   % find the max-min power control
   max_min_power = SINR_candidate * ((IK-SINR_candidate*RhoijMatrix)\b1);
   
   if (a_max - a_min) < diff_compare_threshold && sum(sign(max_min_power)) == n_user && sum(max_min_power) < Ptot && sum(max_min_power) > Ptot - Ptot_margin
       break;
   end
end
iterations = min(index,n_iteration);
%% now, find the SINR of each user, it should be close to SINR_candidate
x = max_min_power;
xnew = x;
SINR_k_maxmin = zeros(1,n_user);
SINR_k_equ = zeros(1,n_user);
for i = 1:n_user
    denom = 1;
    denom_Equ = 1;
    for j = 1:n_user
       if j ~= i
          denom = denom + norm_hi2(i) * rho2_ij(i,j)* xnew(j);
          denom_Equ = denom_Equ + norm_hi2(i) * rho2_ij(i,j)* (Ptot/n_user);
       end
    end
    SINR_k_maxmin(i)  = (norm_hi2(i) * xnew(i))/denom;
    SINR_k_equ(i)     = (norm_hi2(i) * (Ptot/n_user))/denom_Equ;
end
if abs(SINR_candidate - SINR_k_maxmin(1)) > diff_compare_threshold || sum(sign(max_min_power)) ~= n_user || sum(max_min_power) > Ptot+0.001
    error('!!!');
end
Ptot_consumed = sum(max_min_power);
end