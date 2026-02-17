%% Sampling noise
clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

% No of trial for different combinations of values

nitr = 100;
ntrials_list = 50:50:5000;

analyticalStds = zeros(numel(ntrials_list), nitr);
empericalStds  = zeros(numel(ntrials_list), nitr);

for h = 1:numel(ntrials_list)

disp(h)

for itr=1:nitr

orientations     = linspace(0, 180, 10); % linspace(0, 180, 18);
ntrials_per_ori  = ntrials_list(h);
b                = 21;
a                = 0.5; % Does a depend upon b
biasAmp          = 0.5;
sigma_si         = 12;
varGain          = 0.5;
sigma_meta       = 0.2;
Cc               = 0.5; 

% Preallocate arrays
n_theta                  = numel(orientations);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(n_theta, ntrials_per_ori);
theta_resp_all        = zeros(n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s_stim = b + a.*(abs(sind(2*orientations))); %sigma_s_stim = sigma_s_stim';
bias = biasAmp*sind(4*orientations);

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
resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error
resp_err_all_reshaped = resp_err_all(:);

% STD from data
dataErrStd = std(resp_err_all_reshaped);

% Analytical STD
modelParams.orientations        = orientations;
modelParams.b                   = b;
modelParams.a                   = a;
modelParams.biasAmp             = biasAmp;
modelParams.sigma_si            = sigma_si;   
modelParams.varGain             = varGain;
modelParams.Cc                  = Cc;
modelParams.sigma_meta          = sigma_meta;

retData = getAnalyticalSol_MCASANDRE_Contn(modelParams);

analyticalStds(h, itr) = retData.E_sigma_m;
empericalStds(h, itr) = dataErrStd;

end

end


stdEstimationErr = analyticalStds - empericalStds;
x = mean(stdEstimationErr, 2);
y = std(stdEstimationErr, 0, 2);

errorbar(ntrials_list*numel(orientations), ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6);

xlabel("No of trials")
ylabel("(STD data - STD analytical)")
