% close all
% clear all

% Parameters
thetas                  = linspace(0, 350, 18);        % Stimulus orientations (degrees)
ntrials_per_theta       = 200;               
sigma_p_levels          = [1, 2, 3, 4, 5, 6];         % Example: low, medium, high uncertainty (contrast levels)
CV                      = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
sigma_m                 = CV.*sigma_p_levels;           % Meta-uncertainty (log-normal std dev)
Cc                      = 1;                         % Confidence criteria

% Preallocate arrays
n_theta                  = length(thetas);
n_levels                 = length(sigma_p_levels);
total_trials_per_level   = n_theta * ntrials_per_theta;

theta_true_all        = zeros(n_levels, total_trials_per_level);
theta_est_all         = zeros(n_levels, total_trials_per_level);
sigma_p_all           = zeros(n_levels, total_trials_per_level);
sigma_hat_all         = zeros(n_levels, total_trials_per_level);
Vc_all                = zeros(n_levels, total_trials_per_level);
confidence_report_all = zeros(n_levels, total_trials_per_level);
reward_val_all        = zeros(n_levels, total_trials_per_level);
reward_val_all_HC     = zeros(n_levels, total_trials_per_level);
reward_val_all_LC     = zeros(n_levels, total_trials_per_level);
perceptual_error_all  = zeros(n_levels, total_trials_per_level);

% Simulation loop
for level = 1:n_levels
    sigma_p = sigma_p_levels(level);  % Set perceptual noise for this contrast

    for i = 1:n_theta
        theta_true = thetas(i);                        % True orientation
        trials = ntrials_per_theta;

        % Internal estimate (Gaussian noise) - Note: this is not wraped
        % In actual behavioral data this will be wraped
        theta_est = theta_true + sigma_p * randn(trials, 1);
        
        % Subject’s estimate of their uncertainty (meta-uncertainty)
        mu_log = log(sigma_p^2 / sqrt(sigma_m(level)^2 + sigma_p^2));
        sigma_log = sqrt(log(1 + (sigma_m(level)^2 / sigma_p^2)));
        sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
        
        % Confidence variable
        Vc = 1 ./ sigma_hat;
        
        % Store
        range = (i-1)*ntrials_per_theta + 1:i*ntrials_per_theta;

        theta_true_all(level, range)          = theta_true;
        theta_est_all(level, range)           = theta_est;
        sigma_p_all(level, range)             = sigma_p;
        sigma_hat_all(level, range)           = sigma_hat;
        Vc_all(level, range)                  = Vc;
        perceptual_error_all(level, range)    = theta_est - theta_true;
        
        % Confidence report
        confidence_report_all(level, range) = Vc > Cc;
        
        % Update reward values
        err = theta_est - theta_true;
        reward_val_all_HC(level, range) = getReward_HC(err);
        reward_val_all_LC(level, range) = getReward_LC();
        % Reward given considering high and low confidence trials
        reward_val_all(level, range)    = confidence_report_all(level, range).*reward_val_all_HC(level, range) + (1 - confidence_report_all(level, range)).*reward_val_all_LC(level, range); 
    end
end

% Example plot: Confidence vs. error, colored by contrast level
abs_error = abs(perceptual_error_all);

figure; 
subplot(3, 3, 1)
hold on;
cols = lines(length(sigma_p_levels)); % distinct colors

for k = 1:length(sigma_p_levels)
    x = theta_true_all(k, :);
    y = theta_est_all(k, :);

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
    x = theta_true_all(k, :);
    y = Vc_all(k, :);

    % Add jitter if needed
    jitter = (k - 1) * 1.5;
    scatter(x + jitter, y, 1, 'MarkerEdgeColor', cols(k,:), 'DisplayName', ...
        sprintf('\\sigma_p = %.1f', sigma_p_levels(k)));
end

xlabel('True orientation (deg)');
ylabel('Confidence variable');
title('Confidence variable by different contrast level');
legend('Location', 'northwest');
ylim([0, 2000])
grid on;


% Proportion high confidence report
subplot(3, 3, 3)
hold on;
cols = lines(length(sigma_p_levels)); % distinct colors

mean_error = mean(abs_error, 2);
std_error = std(abs_error, 0, 2);

errorbar(linspace(1, numel(sigma_p_levels), numel(sigma_p_levels)), ...
    mean_error, std_error, 'o-', 'LineWidth', 2, 'MarkerSize', 6);

xlabel('Contrast level');
ylabel('Absolute Error');
grid on;
hold off


subplot(3, 3, 4)

LC_mean_errors = zeros(1, length(sigma_p_levels));
HC_mean_errors = zeros(1, length(sigma_p_levels));
LC_std_errors  = zeros(1, length(sigma_p_levels));
HC_std_errors  = zeros(1, length(sigma_p_levels));

for k = 1:length(sigma_p_levels)
    confidenceReports = confidence_report_all(k, :);
    errorReports = abs_error(k, :);

    HC_idx = confidenceReports == 1; % High confidence
    LC_idx = confidenceReports == 0; % Low confidence
    
    errorLC = errorReports(LC_idx);
    errorHC = errorReports(HC_idx);
    
    LC_mean_errors(k) = mean(errorLC);
    HC_mean_errors(k) = mean(errorHC);

    LC_std_errors(k) = std(errorLC);
    HC_std_errors(k) = std(errorHC);
end

hold on
errorbar(linspace(1, numel(sigma_p_levels), numel(sigma_p_levels)), ...
    HC_mean_errors, HC_std_errors, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

errorbar(linspace(1, numel(sigma_p_levels), numel(sigma_p_levels)), ...
    LC_mean_errors, LC_std_errors, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="Low confidence");

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
scatter(sigma_hat_HC(:), errorHC(:), DisplayName="High confidence");
scatter(sigma_hat_LC(:), errorLC(:), DisplayName="Low confidence");
xlabel('Estimated uncertanity');
ylabel('Absolute error (degrees)');
legend();
xlim([0, 100])
hold off

n_bins = 50;
bin_edges = linspace(min(sigma_hat_all(:)), max(sigma_hat_all(:)), n_bins + 1);

% Digitize: assign each sigma_hat_HC value to a bin
[~,~,bin_idx_LC] = histcounts(sigma_hat_LC(:), bin_edges);
[~,~,bin_idx_HC] = histcounts(sigma_hat_HC(:), bin_edges);

% Initialize outputs
mean_error_LC = accumarray(bin_idx_LC, errorLC(:), [n_bins 1], @mean, NaN);
std_error_LC  = accumarray(bin_idx_LC, errorLC(:), [n_bins 1], @std, NaN);

mean_error_HC = accumarray(bin_idx_HC, errorHC(:), [n_bins 1], @mean, NaN);
std_error_HC  = accumarray(bin_idx_HC, errorHC(:), [n_bins 1], @std, NaN);

% Optional: compute bin centers for plotting
bin_centers = 0.5 * (bin_edges(1:end-1) + bin_edges(2:end));

subplot(3, 3, 6)

hold on
errorbar(bin_centers, ...
    mean_error_HC, std_error_HC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");
errorbar(bin_centers, ...
    mean_error_LC, std_error_LC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="Low confidence");
xlabel('Estimated uncertanity');
ylabel('Absolute error (degrees)');
legend();
hold off


% Step 1: Define bins
n_bins = 50;  % or whatever makes sense for your range
edges = linspace(min(sigma_hat_all(:)), max(sigma_hat_all(:)), n_bins+1);
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Step 2: Assign bin indices
[~, ~, bin_idx] = histcounts(sigma_hat_all(:), edges);

% Initialize outputs
mean_reward_LC = accumarray(bin_idx, reward_val_all_LC(:), [n_bins 1], @mean, NaN);
std_reward_LC  = accumarray(bin_idx, reward_val_all_LC(:), [n_bins 1], @std, NaN);

mean_reward_HC = accumarray(bin_idx, reward_val_all_HC(:), [n_bins 1], @mean, NaN);
std_reward_HC  = accumarray(bin_idx, reward_val_all_HC(:), [n_bins 1], @std, NaN);

subplot(3, 3, 7)
hold on
errorbar(bin_centers, ...
    mean_reward_HC, std_reward_HC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");
errorbar(bin_centers, ...
    mean_reward_LC, std_reward_LC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="Low confidence");
xlabel('Estimated uncertanity');
ylabel('Expected reward');
legend();
xlim([0, 100])
hold off

% Expected reward with stimulus uncertainty
mean_reward_LC = mean(reward_val_all_LC, 2);
std_reward_LC = std(reward_val_all_LC, 0, 2);

mean_reward_HC = mean(reward_val_all_HC, 2);
std_reward_HC = std(reward_val_all_HC, 0, 2);

subplot(3, 3, 8)
hold on
errorbar(linspace(1, numel(sigma_p_levels), numel(sigma_p_levels)), ...
    mean_reward_HC, std_reward_HC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");
errorbar(linspace(1, numel(sigma_p_levels), numel(sigma_p_levels)), ...
    mean_reward_LC, std_reward_LC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="Low confidence");
xlabel('Stimulus uncertanity');
ylabel('Expected reward');
legend();
hold off


% Expected reward with stimulus uncertainty
n_bins = 50;  % or whatever makes sense for your range
edges = linspace(min(abs_error(:)), max(abs_error(:)), n_bins+1);
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Step 2: Assign bin indices
[~, ~, bin_idx] = histcounts(abs_error(:), edges);

% Initialize outputs
mean_reward_LC = accumarray(bin_idx, reward_val_all_LC(:), [n_bins 1], @mean, NaN);
std_reward_LC  = accumarray(bin_idx, reward_val_all_LC(:), [n_bins 1], @std, NaN);

mean_reward_HC = accumarray(bin_idx, reward_val_all_HC(:), [n_bins 1], @mean, NaN);
std_reward_HC  = accumarray(bin_idx, reward_val_all_HC(:), [n_bins 1], @std, NaN);

subplot(3, 3, 9)
hold on
errorbar(bin_centers, ...
    mean_reward_HC, std_reward_HC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");
errorbar(bin_centers, ...
    mean_reward_LC, std_reward_LC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="Low confidence");
xlabel('Perceptual error (Absolute)');
ylabel('Expected reward');
legend();
hold off


% Reward function - HC
function [rewardVal] = getReward_HC(error)
    sigma_reward = 1;
    rewardVal = exp( -error.^2 / (2*sigma_reward^2) ) / sqrt(2*pi*sigma_reward^2);
end

function [rewardVal] = getReward_LC()
    rewardVal = 0.05;
end



