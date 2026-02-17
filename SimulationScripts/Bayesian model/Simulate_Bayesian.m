%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

orientations     = 0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 1000;
b                = 0.7;
a                = 0.5; % Does a depend upon b
% biasAmp          = 0.5;
shape            = 2;
scale            = 0.5;
sigma_meta       = 0.2;
Cc               = 1.313;

% Preallocate arrays
n_theta                  = numel(orientations);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(n_theta, ntrials_per_ori);
theta_resp_all        = zeros(n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
% sigma_s_stim = b + a.*(abs(sind(2*orientations))); %sigma_s_stim = sigma_s_stim';
% bias = biasAmp*sind(2*orientations);
% bias = biasAmp*sind(4*orientations);

for i = 1:n_theta

    disp(i)

    theta_true = orientations(i);   % True orientation
    trials = ntrials_per_ori;
    
    % Step1: Define sigma_s_stim
    sigma_s_stim = b + a.*(abs(sind(2*orientations)));
    
    % Step2: Compute sigma_m_stim
    sigma_si = gamrnd(shape, scale, [trials 1]);
    sigma_m_stim = sqrt(sigma_s_stim(i).^2 + sigma_si.^2); %TODO: i
    
    % Step3: Sample estimate of sigma_m from lognormal distribution with spread
    % sigma_meta
    mu_log = log(sigma_m_stim.^2 ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2));
    sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m_stim.^2)));
    sigma_m_stim_hat = lognrnd(mu_log, sigma_log, trials, 1);
    
    % Step4: From the estimate of sigma_m sample orientation
    mu = theta_true;
    sigma = sigma_m_stim_hat;
    theta_measurement = mu + sigma .* randn(trials, 1); 
    theta_measurement = mod(theta_measurement, 180); % Since this is orientation, wrap the angle between 0 and 180
    
    % Step5: MAP estimate - m: theta_measurement, s: theta_true
    % Assuming flat prior p_s_m ∝ p_m_s*p_s
    stimSamplesCnt = 999;
    s = linspace(0, 180, stimSamplesCnt);
    acute_angle_err = angle(theta_measurement, s); % Acute angle error between measurement and true orientation
    % acute_angle_err = theta_measurement - s;
    % Incorporate prior as well
    p_m_s = posterior_prob_dist(acute_angle_err, sigma_m_stim_hat); 
    % exp(-acute_angle_err.^2 ./ (2*sigma_m_stim_hat.^2)) ./ sqrt(2*pi*sigma_m_stim_hat.^2);
    
    [dummy, idx] = max(p_m_s, [], 2);
    theta_estimate = s(idx)';
    
    % raw_err = theta_measurement - s;
    % diff = angle(theta_estimate, theta_measurement);
    
    % Step6: motor noise added on the MAP estimate of measurement
    
    % Step7: Compute confidence
    % confidence variable - noise fluctuations here impacts confidence variable
    % - same as CASANDRE - how is this bayesian model conceptually any
    % different from CASANDRE.
    
    % consider samples of stimuli orientations
    % lognrnd(mu_log, sigma_log) - uniformly sample sigma_m_stim_hat from the
    % lognormal distribution
    % For each sample of sigma_m_stim_hat compute expected error/loss - ds dsigma
    loss = angle(theta_estimate, s).^2;
    
    mu_log = log(sigma_m_stim.^2 ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2)); % For each trial
    sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m_stim.^2)));
    sampleCnt = 500;
    sigma_m_stim_hat_samples  = logninv(linspace(1/sampleCnt, 1 - 1/sampleCnt, sampleCnt), mu_log, sigma_log);
    % sigma_m_stim_hat_samples  = linspace(0.1, 5, sampleCnt); % problamatic
    % sigma_m_stim_hat_samples  = repmat(sigma_m_stim_hat_samples, [ntrials_per_ori, 1]);
    
    confVariables_1 = zeros(ntrials_per_ori, 1);
    confVariables_2 = zeros(ntrials_per_ori, 1);
    
    for nt=1:ntrials_per_ori
        err_this_trial = loss(nt, :);
        sigma_m_samples_this_trial = sigma_m_stim_hat_samples(nt, :)';
        
        s_arr = repmat(s, [sampleCnt, 1]);
        err_arr = repmat(err_this_trial, [sampleCnt, 1]);
        sigma_m_arr = repmat(sigma_m_samples_this_trial, [1, stimSamplesCnt]);
        
        p_ms_ = posterior_prob_dist(err_arr, sigma_m_arr);
        
        % compute confidence
        exp_loss_1 = sum(err_arr.*p_ms_, "all"); % Normalize by ds, dsigma
        
        % Marginalize over s
        ds_ = s(2:end) - s(1:end-1);
        p_ms_1 = ( p_ms_(:,1:end-1) + p_ms_(:,2:end) / 2 ) * ds_';
    
        % Marginalize over dsigma
        dsigma_ = ( sigma_m_samples_this_trial(2:end) - sigma_m_samples_this_trial(1:end-1) ) / 2;
        exp_loss_2 = ( p_ms_1(1:end-1) + p_ms_1(2:end) / 2 )' * dsigma_;
        
        % Expected loss here is same as confidence variable
        confVariables_1(nt) = exp_loss_1;
        confVariables_2(nt) = exp_loss_2;
    end
    
    
    % Set datastructures
    theta_true_all(i, :)         = theta_true;
    theta_resp_all(i, :)         = theta_estimate;
    confidence_report_all(i, :)  = confVariables_2 > Cc;

end

%%

resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all = mod(resp_err_all + 90, 180) - 90;  % Find minimum acute angle error

%% Plot histogram, bias, variance with orientation
%% Histogram (by orientation)
figure 
for i=1:n_theta

subplot(4, 5, i)
histogram(resp_err_all(i, :), Normalization="pdf")
hold on
xline(0, LineStyle="--")
hold off
ylim([0, 1])
xlim([-7, 7])
xlabel("Error (deg)")
ylabel("count")
title(sprintf("Ori %d", orientations(i)))

end

%% Raw err and Confidence (aggregate and by orientation)
figure

resp_err_all_flat = resp_err_all(:);
confidence_report_all_flat = confidence_report_all(:);

resp_err_all_flat_HC = resp_err_all_flat(confidence_report_all_flat == 1);
resp_err_all_flat_LC = resp_err_all_flat(confidence_report_all_flat == 0);

% Average PDF over all the orientation to get analytical solution

subplot(3, 3, 1)
histogram(resp_err_all_flat, Normalization="pdf")
hold on
xlim([-7, 7])
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("All reports, Std: %.2f", std(resp_err_all_flat)))
hold off

subplot(3, 3, 2)
histogram(resp_err_all_flat_HC, Normalization="pdf")
xlim([-7, 7])
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("HC, Std: %.2f", std(resp_err_all_flat_HC)))

subplot(3, 3, 3)
histogram(resp_err_all_flat_LC, Normalization="pdf")
xlim([-7, 7])
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("LC, Std: %.2f", std(resp_err_all_flat_LC)))

propHC = numel(resp_err_all_flat_HC) / numel(resp_err_all_flat);
propLC = numel(resp_err_all_flat_LC) / numel(resp_err_all_flat);

subplot(3, 3, 4)
bar([0 - 0.25, 1 - 0.25], [propHC, propLC], 0.45, DisplayName="data")
xticks([0, 1])
xticklabels({'HC', 'LC'})
legend
ylabel("Proportion (Probability)")

prop_LC_by_ori = arrayfun(@(i) numel(resp_err_all(i, confidence_report_all(i,:) == 0))/size(resp_err_all,2), 1:size(resp_err_all,1))';
prop_HC_by_ori = arrayfun(@(i) numel(resp_err_all(i, confidence_report_all(i,:) == 1))/size(resp_err_all,2), 1:size(resp_err_all,1))';

subplot(3, 3, 5)
hold on
scatter(orientations, prop_HC_by_ori, DisplayName="HC")
scatter(orientations, prop_LC_by_ori, DisplayName="LC")
hold off
xlabel("Orientation (deg)")
ylabel("Proportion (Probability)")
legend

std_err_HC = std(resp_err_all_flat_HC);
std_err_LC = std(resp_err_all_flat_LC);

subplot(3, 3, 6)
% Analytical and actual solution are not exactly equal. Maybe because of
% approximation!
bar([0 - 0.25, 1 - 0.25], [std_err_HC, std_err_LC], 0.45, DisplayName="data")
xticks([0, 1])
xticklabels({'HC', 'LC'})
legend
ylabel("\sigma_m(s)")


% Get LC and HC std per row
std_LC_by_ori = arrayfun(@(i) std(resp_err_all(i, confidence_report_all(i,:) == 0)), 1:size(resp_err_all,1))';
std_HC_by_ori = arrayfun(@(i) std(resp_err_all(i, confidence_report_all(i,:) == 1)), 1:size(resp_err_all,1))';

subplot(3, 3, 7)
hold on
scatter(orientations, std_HC_by_ori, DisplayName='HC')
scatter(orientations, std_LC_by_ori, DisplayName='LC')
hold off
xlabel("Oreintation")
ylabel("\sigma_m(s)")
legend

% Bias 
mean_err   = mean(resp_err_all, 2);
std_err    = std(resp_err_all, [], 2);

subplot(3, 3, 8)
scatter(orientations, mean_err)
hold on
xlabel("orientation (deg)")
ylabel("Mean error")
title("Orientation Bias")
hold off

subplot(3, 3, 9)
bar(orientations, std_err);
hold on
xlabel("orientation (deg)")
ylabel("\sigma_m(s)")
ylim([min(std_err) - 0.1*min(std_err), max(std_err) + 0.1*max(std_err)])
title("Stim dependent measurement noise")
hold off







%%
function [retData] = angle(o1, o2)
err = o1 - o2;
retData = mod(err + 90, 180) - 90;
end

function [retData] = posterior_prob_dist(acute_angle_err, sigma_m)
retData = exp(-acute_angle_err.^2 ./ (2*sigma_m.^2)) ./ sqrt(2*pi*sigma_m.^2);
end
