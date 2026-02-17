%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')
addpath('LL_scripts/')

orientations     = 0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 1000;
b                = 0.9; % Note: different minimum noise level (0.1). Choose b such that average noise level ranges from low to high (relative to internal noise level)
a                = 0.67.*b;   % Does a depend upon b? Yes
% biasAmp          = 0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% scale            = 0.5;
% sigma_meta       = 0.2;
% Cc               = 0.5; 

biasAmp          = 0;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
scale            = 343.3225;
sigma_meta       = 10.83;
Cc               = 0.05; 

% Preallocate arrays
n_theta                  = numel(orientations);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(n_theta, ntrials_per_ori);
theta_resp_all        = zeros(n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
% sigma_s_stim = b + a.*(abs(sind(2*orientations)));
% bias = biasAmp*sind(2*orientations);

sigma_s = [15.5048   24.9028   19.7808   30.5509   34.0230   37.4675]';
sigma_s_stim = b + a.*(abs(sind(2*orientations)));
sigma_s_stim(:) = sigma_s(1);
bias = biasAmp*sind(2*orientations);

for i = 1:n_theta
    theta_true = orientations(i);   % True orientation
    trials = ntrials_per_ori;
    
    % Step1: Using estimate of this uncertainty, the subject estimates
    % the orientation
    % Internal estimate (Gaussian noise) - Note: this is not wraped
    % In actual behavioral data this will be wraped
    % theta_est = theta_true + sigma_p * randn(trials, 1);
    % Find doubly stochastic theta
%     shape = sigma_s_stim(i);
%     gain = gamrnd(shape, scale, [trials 1]);
%     sigma_si_modulated = gain;
%     sigma_m_stim = sqrt(sigma_s_stim(i).^2 + sigma_si_modulated.^2); % For now maybe keep multiplicative noise additive as data seem to show similar trend
%     mean_m_stim = theta_true + bias(i);
    
    % Multiplicative
    shape = sigma_s_stim(i).^2 / scale; % divide by scale so that mean is sigma_s
    gain = gamrnd(shape, scale, [trials 1]);
    sigma_m_stim = sqrt( gain );
    mean_m_stim = theta_true + bias(i);
    
    % TODO: take into account bias?
    % Wrap the angle at the plotting stage. Note: warapping should be
    % performed according to the true angle.
    theta_est = mean_m_stim + sigma_m_stim .* randn(trials, 1);
    % Since this is orientation, wrap the angle between 0 and 180
    theta_est = mod(theta_est, 180); 
    % Are the actual distribution near 0 and 180 captured by this in simulation

    assert(numel(sigma_m_stim) == trials);
    
    % Step1: Subject first gets an estimate of its uncertainty
    % Subject’s estimate of their uncertainty (meta-uncertainty)
    mu_log = log(sigma_m_stim.^2 ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2));
    sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m_stim.^2)));
    sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
    
    % Confidence variable
    Vc = 1 ./ sigma_hat;
    
    % Store
    theta_true_all(i, :)         = theta_true;
    theta_resp_all(i, :)         = theta_est;
    confidence_report_all(i, :)  = Vc > Cc;
end

% Plot the performance curves
resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error

%% Get analytical solution
rvOriErr = -90:0.1:90;

modelParams.sigma_s             = sqrt( mean( sigma_s_stim.^2 ) + std(bias).^2 ); % is this right analytical solition? Maybe not!
modelParams.scale               = scale;
modelParams.Cc                  = Cc;
modelParams.sigma_meta          = sigma_meta;

retData = getEstimationsPDF_cov_reduced(rvOriErr, modelParams);

%% Histogram (by orientation)
figure 
for i=1:n_theta

subplot(4, 5, i)
histogram(resp_err_all(i, :), Normalization="pdf", BinEdges=rvOriErr)
% hold on
% plot(retData.rvOriErrs, retData.analyticalPDF_stim(i, :), LineWidth=1.5);
xline(0, LineStyle="--")
hold off
ylim([0, 1])
xlim([-7, 7])
xlabel("Error (deg)")
ylabel("count")
title(sprintf("Ori %d", orientations(i)))

end

%% HC
figure 
for i=1:n_theta

subplot(4, 5, i)
histogram(resp_err_all(i, confidence_report_all(i, :) == 1), Normalization="pdf")
% hold on
% plot(retData.rvOriErrs, retData.analyticalPDF_stim_HC(i, :), LineWidth=1.5);
xline(0, LineStyle="--")
hold off
ylim([0, 1])
xlim([-7, 7])
xlabel("Error (deg)")
ylabel("count")
title(sprintf("HC: Ori %d", orientations(i)))

end

%% LC
figure 
for i=1:n_theta

subplot(4, 5, i)
histogram(resp_err_all(i, confidence_report_all(i, :) == 0), Normalization="pdf")
% hold on
% plot(retData.rvOriErrs, retData.analyticalPDF_stim_LC(i, :), LineWidth=1.5);
xline(0, LineStyle="--")
hold off
ylim([0, 1])
xlim([-7, 7])
xlabel("Error (deg)")
ylabel("count")
title(sprintf("LC: Ori %d", orientations(i)))

end

%% Raw err and Confidence (aggregate and by orientation)
figure

resp_err_all_flat = resp_err_all(:);
confidence_report_all_flat = confidence_report_all(:);

resp_err_all_flat_HC = resp_err_all_flat(confidence_report_all_flat == 1);
resp_err_all_flat_LC = resp_err_all_flat(confidence_report_all_flat == 0);

subplot(3, 3, 1)
histogram(resp_err_all_flat, Normalization="pdf", BinEdges=rvOriErr)
hold on
plot(retData.rvOriErrs, retData.analyticalPDF, LineWidth=1.5);
xlim([-7, 7])
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("All reports, Std: %.2f, \nAnalytical std: %.2f", std(resp_err_all_flat), retData.E_sigma_m))
hold off

subplot(3, 3, 2)
histogram(resp_err_all_flat_HC, Normalization="pdf", BinEdges=rvOriErr)
hold on
plot(retData.rvOriErrs, retData.analyticalPDF_HC, LineWidth=1.5);
xlim([-7, 7])
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("HC, Std: %.2f, \nAnalytical std: %.2f", std(resp_err_all_flat_HC), retData.E_sigma_m_HC))
hold off

subplot(3, 3, 3)
histogram(resp_err_all_flat_LC, Normalization="pdf", BinEdges=rvOriErr)
hold on
plot(retData.rvOriErrs, retData.analyticalPDF_LC, LineWidth=1.5);
xlim([-7, 7])
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("LC, Std: %.2f, \nAnalytical std: %.2f", std(resp_err_all_flat_LC), retData.E_sigma_m_LC))
hold off

propHC = numel(resp_err_all_flat_HC) / numel(resp_err_all_flat);
propLC = numel(resp_err_all_flat_LC) / numel(resp_err_all_flat);

subplot(3, 3, 4)
bar([0 - 0.25, 1 - 0.25], [propHC, propLC], 0.45, DisplayName="data")
hold on
bar([0 + 0.25, 1 + 0.25], [retData.pHC, retData.pLC], 0.45, DisplayName="Analytical")
hold off
xticks([0, 1])
xticklabels({'HC', 'LC'})
legend
ylabel("Proportion (Probability)")

prop_LC_by_ori = arrayfun(@(i) numel(resp_err_all(i, confidence_report_all(i,:) == 0))/size(resp_err_all,2), 1:size(resp_err_all,1))';
prop_HC_by_ori = arrayfun(@(i) numel(resp_err_all(i, confidence_report_all(i,:) == 1))/size(resp_err_all,2), 1:size(resp_err_all,1))';

subplot(3, 3, 5)
% plot(orientations, retData.pHC_stim, HandleVisibility="off", LineWidth=1.5)
hold on
scatter(orientations, prop_HC_by_ori, DisplayName="HC")
% plot(orientations, retData.pLC_stim, HandleVisibility="off", LineWidth=1.5)
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
hold on
bar([0 + 0.25, 1 + 0.25], [retData.E_sigma_m_HC, retData.E_sigma_m_LC], 0.45, DisplayName="Analytical")
hold off
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
% plot(orientations, retData.E_sigma_m_stim_HC, LineWidth=1.5, HandleVisibility="off")
scatter(orientations, std_LC_by_ori, DisplayName='LC')
% plot(orientations, retData.E_sigma_m_stim_LC, LineWidth=1.5, HandleVisibility="off")
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
% plot(orientations, retData.bias, LineWidth=1.5)
xlabel("orientation (deg)")
ylabel("Mean error")
title("Orientation Bias")
hold off

subplot(3, 3, 9)
bar(orientations, std_err);
hold on
% plot(orientations, retData.E_sigma_m_stim, LineWidth=2)
xlabel("orientation (deg)")
ylabel("\sigma_m(s)")
ylim([min(std_err) - 0.1*min(std_err), max(std_err) + 0.1*max(std_err)])
title("Stim dependent measurement noise")
hold off




