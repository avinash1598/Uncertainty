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
b                = linspace(0.01, 2, 20);
a                = 0.5; % Does a depend upon b
% biasAmp          = 0.5;
shape            = 2;
scale            = 0.5;
sigma_meta       = 0.2;
Cc               = 1;

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
confidence_var_all    = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
% sigma_s_stim = b + a.*(abs(sind(2*orientations))); %sigma_s_stim = sigma_s_stim';
% bias = biasAmp*sind(2*orientations);
% bias = biasAmp*sind(4*orientations);

% Step1: Define sigma_s_stim
sigma_s_stim = b' + a'.*(abs(sind(2*orientations)));      

for l=1:uncertainty_levels
    for i = 1:n_theta
    
        disp(i)
    
        theta_true = orientations(i);   % True orientation
        trials = ntrials_per_ori;
        
        % Step2: Compute sigma_m_stim
        sigma_si = gamrnd(shape, scale, [trials 1]);
        sigma_m_stim = sqrt(sigma_s_stim(l, i).^2 + sigma_si.^2); %TODO: i
        
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
        loss = angle(theta_estimate, s).^2; % take squared error instead
        
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
        theta_true_all(l, i, :)         = theta_true;
        theta_resp_all(l, i, :)         = theta_estimate;
        confidence_report_all(l, i, :)  = confVariables_2 > Cc;
        confidence_var_all(l, i, :)     = confVariables_2;
    
    end
end

%%
% Plot the performance curves
resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error

resp_err_all_reshaped = reshape(resp_err_all, uncertainty_levels, []);
% confidence_report_all_reshaped = reshape(confidence_report_all, uncertainty_levels, []);

cvar_tmp = confidence_var_all > 0.85;
confidence_report_all_reshaped = reshape(cvar_tmp, uncertainty_levels, []);



figure

mean_sigma_s_stim = mean(sigma_s_stim, 2);

x = mean(resp_err_all_reshaped, 2);
y = std(resp_err_all_reshaped, 0, 2);

HC_idx = confidence_report_all_reshaped == 1;
LC_idx = confidence_report_all_reshaped == 0;

resp_HC = resp_err_all_reshaped;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all_reshaped;
resp_LC(~LC_idx) = NaN;

x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');

x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');

x1 = resp_HC(1, :); valid_idx = ~isnan(x1); x1 = x1(valid_idx);
x2 = resp_LC(1, :); valid_idx = ~isnan(x2); x2 = x2(valid_idx);
x3 = resp_HC(uncertainty_levels, :); valid_idx = ~isnan(x3); x3 = x3(valid_idx);
x4 = resp_LC(uncertainty_levels, :); valid_idx = ~isnan(x4); x4 = x4(valid_idx);

subplot(2, 3, 1)
errorbar(mean_sigma_s_stim, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("\sigma_s(s)")
ylabel("Error")

subplot(2, 3, 2)

% Behavioral variability
scatter(mean_sigma_s_stim, y, "filled");
hold on
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim, LineWidth=1.5);
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(mean_sigma_s_stim, y_HC, "filled", DisplayName="High confidence");
hold on
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_HC, LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s_stim, y_LC, "filled", DisplayName="Low confidence");
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_LC, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
legend
hold off




%%
function [retData] = angle(o1, o2)
err = o1 - o2;
retData = mod(err + 90, 180) - 90;
end

function [retData] = posterior_prob_dist(acute_angle_err, sigma_m)
retData = exp(-acute_angle_err.^2 ./ (2*sigma_m.^2)) ./ sqrt(2*pi*sigma_m.^2);
end
