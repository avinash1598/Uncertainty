close all
clear all

% Parameters
thetas = linspace(0, 350, 18);         % Stimulus orientations (degrees)
ntrials_per_theta = 200;               
sigma_p_levels = [1, 3, 5, 7, 9, 11];           % Example: low, medium, high uncertainty (contrast levels)
sigma_m = 6;                           % Meta-uncertainty (log-normal std dev)
Cc = 0.5;                                % Confidence criteria

% Preallocate arrays
n_theta = length(thetas);
n_levels = length(sigma_p_levels);
total_trials = n_theta * ntrials_per_theta * n_levels;

theta_true_all        = zeros(total_trials, 1);
theta_est_all         = zeros(total_trials, 1);
sigma_p_all           = zeros(total_trials, 1);
sigma_hat_all         = zeros(total_trials, 1);
Vc_all                = zeros(total_trials, 1);
contrast_level        = zeros(total_trials, 1);
confidence_report_all = zeros(total_trials, 1);
reward_val_all        = zeros(total_trials, 1);
reward_val_all_HC     = zeros(total_trials, 1);
reward_val_all_LC     = zeros(total_trials, 1);

% Simulation loop
idx = 1;
for level = 1:n_levels
    sigma_p = sigma_p_levels(level);  % Set perceptual noise for this contrast

    for i = 1:n_theta
        theta_true = thetas(i);                        % True orientation
        trials = ntrials_per_theta;

        % Internal estimate (Gaussian noise) - Note: this is not wraped
        % In actual behavioral data this will be wraped
        theta_est = theta_true + sigma_p * randn(trials, 1);
        
        % Subject’s estimate of their uncertainty (meta-uncertainty)
        mu_log = log(sigma_p^2 / sqrt(sigma_m^2 + sigma_p^2));
        sigma_log = sqrt(log(1 + (sigma_m^2 / sigma_p^2)));
        sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
        
        % Confidence variable
        Vc = 1 ./ sigma_hat;
        
        % Store
        range = idx:(idx + trials - 1);
        theta_true_all(range) = theta_true;
        theta_est_all(range)  = theta_est;
        sigma_p_all(range)    = sigma_p;
        sigma_hat_all(range)  = sigma_hat;
        Vc_all(range)         = Vc;
        contrast_level(range) = level;
        
        % Confidence report
        confidence_report_all(range) = Vc > Cc;
        
        % Update reward values
        err = theta_est - theta_true; % How about error beign calculated using sample mean and dispersion - investigate this scenario
        reward_val_all_HC(range) = getReward_HC(err);
        reward_val_all_LC(range) = getReward_LC();
        
        reward_val_all(range) = confidence_report_all(range).*reward_val_all_HC(range) + (1 - confidence_report_all(range)).*reward_val_all_LC(range); 

        idx = idx + trials;
    end
end

% Example plot: Confidence vs. error, colored by contrast level
errors = theta_est_all - theta_true_all;
abs_error = abs(theta_est_all - theta_true_all);

figure; 
subplot(3, 3, 1)
hold on;
cols = lines(length(sigma_p_levels)); % distinct colors

for k = 1:length(sigma_p_levels)
    idx = contrast_level == k;
    x = theta_true_all(idx);
    y = theta_est_all(idx);

    % Add jitter if needed
    jitter = (k - 1) * 1.5;
    scatter(x + jitter, y, 1, 'MarkerEdgeColor', cols(k,:), 'DisplayName', ...
        sprintf('\\sigma_p = %.1f', sigma_p_levels(k)));
end

xlabel('True orientation (deg)');
ylabel('Estimated orientation (deg)');
title('Estimated vs. True Orientation by Contrast Level');
legend('Location', 'northwest');
grid on;


subplot(3, 3, 2)
hold on;
cols = lines(length(sigma_p_levels)); % distinct colors

for k = 1:length(sigma_p_levels)
    idx = contrast_level == k;
    x = theta_true_all(idx);
    y = Vc_all(idx);

    % Add jitter if needed
    jitter = (k - 1) * 1.5;
    scatter(x + jitter, y, 1, 'MarkerEdgeColor', cols(k,:), 'DisplayName', ...
        sprintf('\\sigma_p = %.1f', sigma_p_levels(k)));
end

xlabel('True orientation (deg)');
ylabel('Confidence variable');
title('Confidence variable by different contrast level');
legend('Location', 'northwest');
grid on;


% Proportion high confidence report
subplot(3, 3, 3)
hold on;
cols = lines(length(sigma_p_levels)); % distinct colors

for k = 1:length(sigma_p_levels)
    idx = contrast_level == k;
    x = abs_error(idx);
    mean_error = mean(x);
    
    % Add jitter if needed
    scatter(k, mean_error, 'Color', cols(k,:), 'DisplayName', ...
        sprintf('\\sigma_p = %.1f', sigma_p_levels(k)));
end

xlabel('Contrast level');
ylabel('Mean Error');
legend('Location', 'northwest');
grid on;
hold off


subplot(3, 3, 4)

LC_errors = zeros(1, length(sigma_p_levels));
HC_errors = zeros(1, length(sigma_p_levels));

for k = 1:length(sigma_p_levels)
    idx = contrast_level == k;
    
    confidenceReports = confidence_report_all(idx);
    errorReports = abs_error(idx);

    HC_idx = confidenceReports == 1; % High confidence
    LC_idx = confidenceReports == 0; % Low confidence
    
    errorLC = errorReports(LC_idx);
    errorHC = errorReports(HC_idx);
    
    LC_errors(k) = mean(errorLC);
    HC_errors(k) = mean(errorHC);
end

hold on
plot(HC_errors, DisplayName="High confidence");
plot(LC_errors, DisplayName="Low confidence");
xlabel('Stimulus uncertanity');
ylabel('Absolute error (degrees)');
legend();
hold off


subplot(3, 3, 5)

HC_idx = confidence_report_all == 1; % High confidence
LC_idx = confidence_report_all == 0; % Low confidence

errorLC = abs_error(LC_idx);
errorHC = abs_error(HC_idx);

sigma_hat_LC = sigma_hat_all(LC_idx);
sigma_hat_HC = sigma_hat_all(HC_idx);

hold on
scatter(sigma_hat_HC, errorHC, DisplayName="High confidence");
scatter(sigma_hat_LC, errorLC, DisplayName="Low confidence");
xlabel('Estimated uncertanity');
ylabel('Absolute error (degrees)');
legend();
hold off


subplot(3, 3, 6)
hold on
scatter(sigma_hat_all, reward_val_all_HC)
scatter(sigma_hat_all, reward_val_all_LC)
xlabel("Estimated uncertainty")
ylabel("Reward")
xlim([0, 50])
hold off            


subplot(3, 3, 7)

HC_idx = confidence_report_all == 1; % High confidence
LC_idx = confidence_report_all == 0; % Low confidence

rewardHC = reward_val_all_HC;
rewardLC = reward_val_all_LC;

hold on
scatter(sigma_hat_all, rewardHC)
scatter(sigma_hat_all, rewardLC)
xlabel("Estimated uncertainty")
ylabel("Reward")
xlim([0, 50])
hold off  


% Step 1: Define bins
n_bins = 50;  % or whatever makes sense for your range
edges = linspace(min(sigma_hat_all), max(sigma_hat_all), n_bins+1);
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Step 2: Assign bin indices
[~, ~, bin_idx] = histcounts(sigma_hat_all, edges);

% Step 3: Compute average reward per bin
avg_reward_per_bin_HC = zeros(n_bins, 1);
avg_reward_per_bin_LC = zeros(n_bins, 1);

for i = 1:n_bins
    this_bin_reward_HC = reward_val_all_HC(bin_idx == i);
    this_bin_reward_LC = reward_val_all_LC(bin_idx == i);
    
    if ~isempty(this_bin_reward_HC)
        avg_reward_per_bin_HC(i) = mean(this_bin_reward_HC);
    else
        avg_reward_per_bin_HC(i) = NaN;  % in case a bin has no values
    end
    
    if ~isempty(avg_reward_per_bin_LC)
        avg_reward_per_bin_LC(i) = mean(this_bin_reward_LC);
    else
        avg_reward_per_bin_LC(i) = NaN;
    end
end

subplot(3, 3, 9)
hold on
plot(bin_centers, avg_reward_per_bin_HC, DisplayName="HC")
plot(bin_centers, avg_reward_per_bin_LC, DisplayName="LC")
xlabel("Estimated uncertainty")
ylabel("Expected reward")
xlim([0, 100])
legend()
hold off


subplot(3, 3, 8)
hold on
scatter(sigma_hat_all, abs_error)
xlabel("Estimated uncertainty")
ylabel("Error")
xlim([0, 50])
hold off

% TODO: add reward scheme
% Do we need to consider orientation dependence

% figure
% scatter(abs_error, reward_val_all)
% xlabel("Error")
% ylabel("Reward")
% 

% Reward function - HC
function [rewardVal] = getReward_HC(error)
    sigma_reward = 5;
    rewardVal = exp( -error.^2 / (2*sigma_reward^2) ) / sqrt(2*pi*sigma_reward^2);
end

function [rewardVal] = getReward_LC()
    rewardVal = 0.01;
end



