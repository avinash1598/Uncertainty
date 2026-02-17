%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

orientations     = 0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 5000;
b                = linspace(0.01, 2, 20);
a                = 0.2;   % Does a depend upon b?
biasAmp          = 0.5;   % Does bias depend upon uncertainty level?
sigma_si         = 1;
varGain          = 0.5;
sigma_meta       = 0.2;
Cc               = 0.5; 

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s_stim = b' + a'.*(abs(sind(2*orientations))); %sigma_s_stim = sigma_s_stim';
bias = biasAmp*sind(2*orientations); % Does bias depend upon uncertainty level?
% bias = biasAmp*sind(4*orientations);

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
        gain = gamrnd(1./varGain, varGain, [trials 1]);
        sigma_si_modulated = gain*sigma_si;
        sigma_m_stim = sqrt(sigma_s_stim(l, i).^2 + sigma_si_modulated.^2);
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

%% Get analytical solution
anlytcl_sigma_m_stim = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_LC = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    modelParams.orientations        = orientations;
    modelParams.b                   = b(i);
    modelParams.a                   = a;
    modelParams.biasAmp             = biasAmp;
    modelParams.sigma_si            = sigma_si;   
    modelParams.varGain             = varGain;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    
    retData = getAnalyticalSol_MCASANDRE_Contn(modelParams);

    anlytcl_sigma_m_stim(i)    = retData.E_sigma_m;
    anlytcl_sigma_m_stim_HC(i) = retData.E_sigma_m_HC;
    anlytcl_sigma_m_stim_LC(i) = retData.E_sigma_m_LC;
end

%% Plot results

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

