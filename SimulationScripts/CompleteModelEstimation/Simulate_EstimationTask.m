%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')
addpath('./LL_scripts/')

orientations     = 0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 1000;
b                = 0.1;
a                = 0.67*b; % Does a depend upon b
biasAmp          = 0.5;
scale            = 0.5;
shape            = 2;   % Change value to set different mean level noise
sigma_meta       = 0.2;
Cc               = 0.5; %0.5

% Preallocate arrays
n_theta                  = numel(orientations);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(n_theta, ntrials_per_ori);
theta_resp_all        = zeros(n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s_stim = b + a.*(abs(sind(2*orientations))); %sigma_s_stim = sigma_s_stim';
bias = biasAmp*sind(2*orientations);
% bias = biasAmp*sind(4*orientations);

for i = 1:n_theta
    theta_true = orientations(i);   % True orientation
    trials = ntrials_per_ori;
    
    % Step1: Using estimate of this uncertainty, the subject estimates
    % the orientation
    % Internal estimate (Gaussian noise) - Note: this is not wraped
    % In actual behavioral data this will be wraped
    % theta_est = theta_true + sigma_p * randn(trials, 1);
    % Find doubly stochastic theta
    gain = gamrnd(shape, scale, [trials 1]);
    sigma_si_modulated = gain;
    sigma_m_stim = sqrt(sigma_s_stim(i).^2 + sigma_si_modulated.^2);
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

% Find minimum acute angle error
resp_err_all = mod(resp_err_all + 90, 180) - 90;

%% Get analytical solution
rvOriErr = -90:0.1:90;
% modelParams.orientations        = orientations;
modelParams.b                   = b;
modelParams.a                   = a;
modelParams.biasAmp             = biasAmp;
modelParams.shape               = shape;   
modelParams.scale               = scale;
modelParams.Cc                  = Cc;
modelParams.sigma_meta          = sigma_meta;

retData = getAnalyticalSol_EstimationTask(orientations, rvOriErr, modelParams);

%% statistical test
A  = resp_err_all; resp_err_by_level = reshape(A, size(A, 1), []); 
A  = confidence_report_all; confidence_report_by_level   = reshape(A, size(A, 1), []);

err_HC_1 = resp_err_by_level(1, confidence_report_by_level(1, :) == 1);
err_HC_2 = resp_err_by_level(6, confidence_report_by_level(6, :) == 1);
err_LC_1 = resp_err_by_level(1, confidence_report_by_level(1, :) == 0);
err_LC_2 = resp_err_by_level(6, confidence_report_by_level(6, :) == 0);

d1 = std(err_LC_1) - std(err_HC_1);
d2 = std(err_LC_2) - std(err_HC_2);
D_obs = d2 - d1;

nBoot = 10000;
D_boot = zeros(nBoot,1);

% Do this with simulated data

for i = 1:nBoot
    sHC1 = std(err_HC_1(randi(numel(err_HC_1), [numel(err_HC_1), 1])));
    sLC1 = std(err_LC_1(randi(numel(err_LC_1), [numel(err_LC_1), 1])));
    sHC2 = std(err_HC_2(randi(numel(err_HC_2), [numel(err_HC_2), 1])));
    sLC2 = std(err_LC_2(randi(numel(err_LC_2), [numel(err_LC_2), 1])));
    
    D_boot(i) = (sLC2 - sHC2) - (sLC1 - sHC1);
end

p = mean(D_boot >= D_obs);

% 95% percentile CI
ci = prctile(D_boot, [2.5 97.5]);

% one-sided p estimate (probability D_boot <= 0)
p_one_sided = mean(D_boot <= 0); % small means evidence D>0

fprintf('D_obs=%.4f, CI=[%.4f, %.4f], p_one_sided (D>0) ~= %.4f\n', D_obs, ci(1), ci(2), p_one_sided);

figure
hold on
histogram(D_boot)
xline(D_obs, LineStyle="--")
xline(ci(1), LineStyle="-")
xline(ci(2), LineStyle="-")
hold off

retData = doBootstrapTest(err_HC_1, err_LC_1, err_HC_2, err_LC_2, false);
retData

%% GLME

n_uncertainty_levels  = uncertainty_levels;
uncertainty_level_all = repmat((1:n_uncertainty_levels)', 1, size(resp_err_by_level, 2));
error                 = abs(resp_err_by_level(:)) ; %+ 1e-12;
uncertaintyLevel      = uncertainty_level_all(:);
confidence            = confidence_report_by_level(:);
T                     = table(error, uncertaintyLevel, confidence);
T.confidence          = categorical(T.confidence);  % HC vs LC

% This analysis might be okay if i correct for bias
glme = fitglme(T, ...
    'error ~ uncertaintyLevel * confidence + (1|uncertaintyLevel)', ...
    'Distribution','Gamma', ...
    'Link','log');

disp(glme)


%% Histogram (by orientation)
figure 
for i=1:n_theta

subplot(4, 5, i)
histogram(resp_err_all(i, :), Normalization="pdf", BinEdges=rvOriErr)
hold on
plot(retData.rvOriErrs, retData.analyticalPDF_stim(i, :), LineWidth=1.5);
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
hold on
plot(retData.rvOriErrs, retData.analyticalPDF_stim_HC(i, :), LineWidth=1.5);
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
hold on
plot(retData.rvOriErrs, retData.analyticalPDF_stim_LC(i, :), LineWidth=1.5);
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
plot(orientations, retData.pHC_stim, HandleVisibility="off", LineWidth=1.5)
hold on
scatter(orientations, prop_HC_by_ori, DisplayName="HC")
plot(orientations, retData.pLC_stim, HandleVisibility="off", LineWidth=1.5)
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
plot(orientations, retData.E_sigma_m_stim_HC, LineWidth=1.5, HandleVisibility="off")
scatter(orientations, std_LC_by_ori, DisplayName='LC')
plot(orientations, retData.E_sigma_m_stim_LC, LineWidth=1.5, HandleVisibility="off")
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
plot(orientations, retData.bias, LineWidth=1.5)
xlabel("orientation (deg)")
ylabel("Mean error")
title("Orientation Bias")
hold off

subplot(3, 3, 9)
bar(orientations, std_err);
hold on
plot(orientations, retData.E_sigma_m_stim, LineWidth=2)
xlabel("orientation (deg)")
ylabel("\sigma_m(s)")
ylim([min(std_err) - 0.1*min(std_err), max(std_err) + 0.1*max(std_err)])
title("Stim dependent measurement noise")
hold off




