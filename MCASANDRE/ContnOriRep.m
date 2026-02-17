close all
clear all

% Simulation parameters
n_trials = 1000;
theta_mean = 350;                % Stimulus range in degrees
sigma_p = 5;                     % Perceptual uncertainty (degrees) - add many levels based on contrast
sigma_m = 2;                     % Meta-uncertainty (log-normal std dev)

% Draw from normal distribution
theta_samples = theta_mean + sigma_p * randn(n_trials, 1);

% Wrap to circular range [0, 360)
theta_wrapped = mod(theta_samples, 360);

% Subject's noisy estimate of their perceptual uncertainty (meta-uncertainty)
% Find log normal distribution parameters
param_lognd_mu = log(sigma_p^2 / sqrt(sigma_m^2 + sigma_p^2));
param_lognd_sigma = sqrt(log(1 + (sigma_m^2 / sigma_p^2)));
sigma_hat = lognrnd(param_lognd_mu, param_lognd_sigma, n_trials, 1);

% Compute confidence variable (inverse of estimated uncertainty)
Vc = 1 ./ sigma_hat;

% Plot
figure;

subplot(1, 2, 1)
histogram(theta_samples, 30);
xlabel('Orientation (degrees)');
ylabel('Count');
title('Wrapped Normal Samples: Centered at 10°');

subplot(1, 2, 2)
histogram(sigma_hat, 50);
xlabel('Confidence variable (Vc)');
ylabel('Frequency');
title('Distribution of confidence variable (Vc)');
