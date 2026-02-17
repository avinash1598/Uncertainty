%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('LL_scripts/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

orientations     = 0:15:175; %linspace(0, 180, 18); %0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 25; %250;
b                = linspace(1, 2.2, 6); % linspace(1, 2.2, 8); Note: different minimum noise level (0.1). Choose b such that average noise level ranges from low to high (relative to internal noise level)
a                = 0.67.*b; %0.67.*b;   % Does a depend upon b? Yes
biasAmp          = 0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
scale            = 0.5; %0.5;
sigma_meta       = 0.6;
Cc               = 0.5; 

% In actual data correct for bias

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s_stim = b' + a'*(abs(sind(2*orientations)));
bias = biasAmp*sind(2*orientations); 

for l=1:uncertainty_levels
    for i = 1:n_theta
        theta_true = orientations(i);   % True orientation
        trials = ntrials_per_ori;
        
        % Step1: Using estimate of this uncertainty, the subject estimates
        % the orientation
        % Internal estimate (Gaussian noise) - Note: this is not wraped
        % In actual behavioral data this will be wraped
        % theta_est = theta_true + sigma_p * randn(trials, 1);
        % Find doubly stochastic theta
%         shape = sigma_s_stim(l, i);
%         gain = gamrnd(shape, scale, [trials 1]);
%         sigma_si_modulated = gain;
%         sigma_m_stim = sqrt(sigma_s_stim(l, i).^2 + sigma_si_modulated.^2); % For now maybe keep multiplicative noise additive as data seem to show similar trend
%         mean_m_stim = theta_true + bias(i);
        
        % Multiplicative
        shape = sigma_s_stim(l, i).^2 / scale; % divide by scale so that mean is sigma_s
        gain = gamrnd(shape, scale, [trials 1]);
        sigma_m_stim = sqrt( gain );
        mean_m_stim = theta_true + bias(i);
        
        % Poisson noise
        % lambda = sigma_s_stim(l, i);
        % sigma_m_stim = poissrnd(lambda, [trials 1]);
        % mean_m_stim = theta_true + bias(i);
        
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
        theta_true_all(l, i, :)         = theta_true;
        theta_resp_all(l, i, :)         = theta_est;
        confidence_report_all(l, i, :)  = Vc > Cc;
    end
end

% Plot the performance curves
resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error

resp_err_all_reshaped = reshape(resp_err_all, uncertainty_levels, []);
confidence_report_all_reshaped = reshape(confidence_report_all, uncertainty_levels, []);

% Save model data
data.stimOri     = theta_true_all;
data.reportedOri = theta_resp_all;
data.err         = resp_err_all;
data.confReport  = confidence_report_all;

data.params.sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
data.params.b                     = b;
data.params.a                     = a;
data.params.biasAmp               = biasAmp;
data.params.scale                 = scale;
data.params.sigma_meta            = sigma_meta;
data.params.Cc                    = Cc;

save('modelContOriData_cov.mat', "data")


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

retData = doBootstrapTest(abs(err_HC_1), abs(err_LC_1), abs(err_HC_2), abs(err_LC_2), true);
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

%% Get analytical solution
anlytcl_sigma_m_stim = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_LC = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    rvOriErr = -90:0.1:90;
    
    modelParams.b                   = b(i);
    modelParams.a                   = a(i);
    modelParams.biasAmp             = biasAmp;
    modelParams.scale               = scale;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    modelParams.guessRate           = 0;
    
    retData = getEstimationsPDF_cov(orientations, rvOriErr, modelParams);
    
    anlytcl_sigma_m_stim(i)    = retData.E_sigma_m;
    anlytcl_sigma_m_stim_HC(i) = retData.E_sigma_m_HC;
    anlytcl_sigma_m_stim_LC(i) = retData.E_sigma_m_LC;
end

%% Plot results

figure

mean_sigma_s_stim = mean(sigma_s_stim, 2);

x = mean(resp_err_all_reshaped, 2);
y = std(resp_err_all_reshaped, 0, 2);
% y = mad(resp_err_all_reshaped, 1, 2);

HC_idx = confidence_report_all_reshaped == 1;
LC_idx = confidence_report_all_reshaped == 0;

resp_HC = resp_err_all_reshaped;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all_reshaped;
resp_LC(~LC_idx) = NaN;

x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');
% y_HC = mad(resp_HC, 1, 2);

x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');
% y_LC = mad(resp_LC, 1, 2);

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
plot(mean_sigma_s_stim, anlytcl_sigma_m_stim, LineWidth=1.5);
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(mean_sigma_s_stim, y_HC, "filled", DisplayName="High confidence");
hold on
plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_HC, LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s_stim, y_LC, "filled", DisplayName="Low confidence");
plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_LC, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
legend
hold off

