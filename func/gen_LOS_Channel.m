function [channel_unit_norm,channel_gain,path_loss_dB] = gen_LOS_Channel(M,K,f0_gig,R_max,R_min,...
                    alpha_LOS,theta_min, theta_max, min_spacing_distance_user, spacing_array, shadowing, flag_correlation, separation_angle)
% Initialization:
n_user = K;             % # single-antenna user
n_bs = M;               % # bs antennas
gig = 10^9;             % giga = 10^9
cwave = 0.3 * gig;      % C: the speed of light
freq = f0_gig;          % the frequency in Ghz
w = 2*pi*freq;          % w = 2*pi*f
k = w/cwave;            % the wave number
lambda = cwave/freq;    % the wavelenght
% the minimum allowable spacing between the users in x dimension
min_x_spacing = max(lambda,min_spacing_distance_user);
% the minimum allowable spacing between the users in y dimension
min_y_spacing = max(lambda,min_spacing_distance_user);
% the coordination of the antennas in the linear array
x_ant_coordinate = spacing_array * (lambda)*(1-n_bs/2:n_bs/2);
%% Distributing the users uniformly in a sector
% distributed uniformly in R^2, and uniformly in phi
flag_user_spacing = 1;  % flag is true as long as the users do not meet the min spacing
% The while block for distributing the users
while flag_user_spacing == 1
    % uniform distribution in angles
    % phi \in [\phi_min,\phi_max]
    phi_user_deg = theta_min + (theta_max-theta_min).*rand(1,n_user);
    if flag_correlation == 1
        phi_user_deg(2) = 90;
        phi_user_deg(1) = phi_user_deg(2) + separation_angle;
        phi_user_deg(3) = theta_min;
    end
    % uniform distribution in R^2
    % generate a uniform random variable between [0,1]
    % then use "a" to generate R_user
    a = (R_min/R_max)^2;
    add_val = a + (1-a)*rand(1,n_user);
    R_user = R_max * sqrt(add_val);
    % Sort the R to further check the min spacing
    [~,index_R] = sort(R_user,'ascend');
    % The location is found by:
    x_user = R_user .* cosd(phi_user_deg);
    y_user = R_user .* sind(phi_user_deg);
    % just sort the x axis to check the minimum spacing distance
    x_user_sort_R = x_user;
    y_user_sort_R = y_user;
    [x_user,index_x_user] = sort(x_user);
    y_user = y_user(index_x_user);
    %% Check whether two users are very close to each other:
    % when both the dimensions of the users are closer than the min
    % spacing, then the users are really really close!
    % we avoid these scenarios! by repeating the generating the users
    x_user_difference_flag =  abs(diff(x_user)) <= min_x_spacing;
    y_user_difference_flag =  abs(diff(y_user)) <= min_y_spacing;
    if sum(x_user_difference_flag.*y_user_difference_flag) == 0
    	flag_user_spacing = 0;
    else
    	flag_user_spacing = 1;
    end
    % Just sort the users based on the R, closer users be first indexes
    [~,index_R] = sort(R_user,'ascend');
    x_user = x_user_sort_R(index_R);
    y_user = y_user_sort_R(index_R);
end
%% We use ETSI TR 138 901 V14.0.0 for the large-scale fading
% and the path loss, and cell configuration
% UMi- Street Canyon is for LOS
%% Build the LOS and NLOS components
channel_nonnormalized   = zeros(n_bs,n_user);
channel_unit_norm       = zeros(n_bs,n_user);
channel_gain            = zeros(1,n_user);
path_loss_dB            = zeros(1,n_user);
%% For LOS channel we use "UMi- Street Canyon" from 3GPP
% the followings assumptions are based on ETSI TR 138 901 V14.0.0 (2017-05)
% UMi-Street Canyon
% ISD = 200
% Min BS-UT = 10m (2D distance)
h_UT = 1.5;
h_BS = 10;
h_E = 1; % for UMi
heffective_BS = h_BS - h_E;
heffective_UT = h_UT - h_E;
% f0, cwave
d_BP = 4*heffective_BS*heffective_UT*(f0_gig/cwave);
%% Main loop for the LOS path loss
if alpha_LOS == 1
    rho_SF = 4; % look at table Table 7.4.1-1: Pathloss models
    % if there is shadowing, use the following code
    for l = 1:n_user
        if shadowing == 1
            lognormal_attenuation = rho_SF*randn(1);
        else
            lognormal_attenuation = 0;
        end
        % see Figure7.4.1-1 of ETSI TR 138 901 V14.0.0
        % for the followings
        % the center of the array locates in the origin
        % thus, d_2D = sqrt(xuser^2 + yuser^2)
        d_2D = sqrt((x_user(l))^2 + (y_user(l))^2);
        d_3D = sqrt((h_UT-h_BS)^2 + d_2D^2 );
        if d_2D <= d_BP
            path_loss_dB(l) = 32.4 +21*log10(d_3D) + 20*log10(f0_gig/gig) + lognormal_attenuation;
        else
            path_loss_dB(l) = 32.4 +40*log10(d_3D) + 20*log10(f0_gig/gig) -9.5*log10(d_BP^2 + (h_BS-h_UT)^2) +lognormal_attenuation;
        end
        % convert dB to decimal
        path_loss_dec = 10^(-path_loss_dB(l)/10);
        for i = 1:n_bs
            R_calculated_LOS           = sqrt((h_UT - h_BS)^2 +(y_user(l))^2 + ( (x_ant_coordinate(i)-x_user(l))^2 ) );
            channel_nonnormalized(i,l) = sqrt(path_loss_dec) * exp(-1j*k*R_calculated_LOS);
        end
        channel_gain(l)         = norm(channel_nonnormalized(:,l));
        channel_unit_norm(:,l)  = channel_nonnormalized(:,l)/channel_gain(l);
    end
else
    % for NLOS
    rho_SF_NLOS = 8.2;  % look at table Table 7.4.1-1: Pathloss models
    for l = 1:n_user
        if shadowing == 1
            lognormal_attenuation = rho_SF_NLOS*randn(1);
        else
            lognormal_attenuation = 0;
        end
        % see Figure7.4.1-1 of ETSI TR 138 901 V14.0.0
        % for the followings
        d_2D = sqrt((x_user(l))^2 + (y_user(l))^2);
        d_3D = sqrt((h_UT-h_BS)^2 + d_2D^2 );
        path_loss_dB(l) = 32.4 + 20*log10(f0_gig/gig) + 31.9*log10(d_3D) + lognormal_attenuation;
        % convert dB to decimal
        path_loss_dec = 10^(-path_loss_dB(l)/10);
        for i = 1:n_bs
            R_calculated_LOS           = sqrt((h_UT - h_BS)^2 +(x_user(l))^2 + ( (x_ant_coordinate(i)-y_user(l))^2 ) );
            channel_nonnormalized(i,l) = sqrt(path_loss_dec) * exp(-1j*k*R_calculated_LOS);
        end
        channel_gain(l)         = norm(channel_nonnormalized(:,l));
        channel_unit_norm(:,l)  = channel_nonnormalized(:,l)/channel_gain(l);
    end
end
end