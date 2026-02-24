%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')

orientations     = linspace(0, 179, 10); %0:10:180; % 
ntrials_per_ori  = 25; %1000;
b                = linspace(0.1, 1.5, 6); % 1.2 % Choose b such that average noise level ranges from low to high (relative to internal noise level)
a                = 0.67.*b; %0.67   % Does a depend upon b? Yes
biasAmp          = 0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
shape            = 2;
scale            = 0.5;
sigma_meta       = 0.2;
Cc               = 0.7; 
guessRate        = 0.1;

% biasAmp          = 0; %0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% shape            = 0.0848; %2;
% scale            = 338.1997; %0.5;
% sigma_meta       = 41.5023; %0.2;
% Cc               = 0.1097; %0.5; 

% biasAmp          = 0;
% shape            = 2; %2;
% scale            = 500; %0.5;
% sigma_meta       = 20; %0.2;
% Cc               = 0.05; %0.5; 

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(b);
n_uncertainty_levels     = numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
% sigma_s = [8.0253, 13.1054, 13.4657, 25.9641, 30.2290, 36.0350]';
% sigma_s = [7.7131, 11.268, 16.754, 25.546, 41.216, 45.754 ]';
% sigma_s_stim = repmat(sigma_s, [1 n_theta]);
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
        gain = gamrnd(shape, scale, [trials 1]);
        sigma_si_modulated = gain;
        sigma_m_stim = sqrt(sigma_s_stim(l, i).^2 + sigma_si_modulated);
        mean_m_stim = theta_true + bias(i);
        
        % TODO: take into account bias?
        % Wrap the angle at the plotting stage. Note: warapping should be
        % performed according to the true angle.
        theta_est = mean_m_stim + sigma_m_stim .* randn(trials, 1);
        theta_est = mod(theta_est, 180); 
        
        % Simulat guess rate
        guess_tl_idx = randi([1 trials], floor( trials*guessRate ), 1);
        guessOris = 180*rand(numel(guess_tl_idx), 1);
        theta_est(guess_tl_idx) = guessOris;
        
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
data.stimOri                = theta_true_all;
data.reportedOri            = theta_resp_all;
data.resp_err_all           = resp_err_all;
data.confidence_report_all  = confidence_report_all;
% data.err         = resp_err_all;
% data.confReport  = confidence_report_all;

data.params.sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
data.params.b                     = b;
data.params.a                     = a;
data.params.biasAmp               = biasAmp;
data.params.shape                 = shape;
data.params.scale                 = scale;
data.params.sigma_meta            = sigma_meta;
data.params.Cc                    = Cc;
data.params.guessRate             = guessRate;

%% Optimize
errBins = -90:0.5:90;

optParams.nStarts = 10;
optParams.hyperParamC1 = 0;
optParams.randomGuessModel = true;

result = Optimize(data, errBins, "ind", [], optParams, "full");

%%
[~, idx] = min(result.f);

opt_param_sigma_s         = result.x(idx, 1:n_uncertainty_levels);
opt_param_shape           = result.x(idx, n_uncertainty_levels + 1 - 0);
opt_param_scale           = result.x(idx, n_uncertainty_levels + 2-0);
opt_param_sigma_meta      = result.x(idx, n_uncertainty_levels + 3-0);
opt_param_Cc              = result.x(idx, n_uncertainty_levels + 4-0);
opt_param_guessrate       = result.x(idx, n_uncertainty_levels + 5-0);
opt_param_sigma_ori_scale = result.x(idx, n_uncertainty_levels + 6-0);
opt_param_bias            = result.x(idx, n_uncertainty_levels + 7-0);

% opt_param_sigma_s         = result.x(idx, 1:n_uncertainty_levels);
% opt_param_scale           = result.x(idx ,n_uncertainty_levels + 1);
% opt_param_sigma_meta      = result.x(idx, n_uncertainty_levels + 2);
% opt_param_Cc              = result.x(idx, n_uncertainty_levels + 3);
% opt_param_guessrate       = result.x(idx, n_uncertainty_levels + 4);
% opt_param_sigma_ori_scale = result.x(idx, n_uncertainty_levels + 5);
% opt_param_bias            = result.x(idx, n_uncertainty_levels + 6);

gt_sigma_s          = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 );
gt_shape            = shape;
gt_scale            = scale;
gt_sigma_meta       = sigma_meta;
gt_Cc               = Cc;
gt_guessrate        = guessRate;
gt_sigma_ori_scale  = mean( a/b );
gt_bias             = biasAmp;

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("GT: %.4f, Fit: %.4f \n", b(i), opt_param_sigma_s(i))
end

fprintf("GT: %.4f, Fit: %.4f \n", gt_shape, opt_param_shape)
fprintf("GT: %.4f, Fit: %.4f \n", gt_scale, opt_param_scale)
fprintf("GT: %.4f, Fit: %.4f \n", gt_sigma_meta, opt_param_sigma_meta)
fprintf("GT: %.4f, Fit: %.4f \n", gt_Cc, opt_param_Cc)
fprintf("GT: %.4f, Fit: %.4f \n", gt_guessrate, opt_param_guessrate)
fprintf("GT: %.4f, Fit: %.4f \n", gt_sigma_ori_scale, opt_param_sigma_ori_scale)
fprintf("GT: %.4f, Fit: %.4f \n", gt_bias, opt_param_bias)

%% Get analytical solution
anlytcl_sigma_m_stim = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_LC = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim_LC = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    rvOriErr = errBins;
    
    modelParams.b                   = opt_param_sigma_s(i);
    modelParams.a                   = gt_sigma_ori_scale*opt_param_sigma_s(i);
    modelParams.biasAmp             = opt_param_bias;
    modelParams.shape               = opt_param_shape;
    modelParams.scale               = opt_param_scale;
    modelParams.Cc                  = opt_param_Cc;
    modelParams.sigma_meta          = opt_param_sigma_meta;
    modelParams.guessRate           = opt_param_guessrate;
    
    % retData = getEstimatesPDFs(1:10:180, rvOriErr, modelParams);
    retData = getEstimationsPDF_cov(1:10:180, rvOriErr, modelParams);
    
    anlytcl_sigma_m_stim(i)    = retData.E_sigma_m;
    anlytcl_sigma_m_stim_HC(i) = retData.E_sigma_m_HC;
    anlytcl_sigma_m_stim_LC(i) = retData.E_sigma_m_LC;

    anlytcl_mad_m_stim(i)      = retData.mad_m;
    anlytcl_mad_m_stim_HC(i)   = retData.mad_m_HC;
    anlytcl_mad_m_stim_LC(i)   = retData.mad_m_LC;

end

%%
% PDFs by uncertainty
figure

n_uncertainty_levels = numel(b);

for i=1:n_uncertainty_levels
    
    modelParams.b                   = opt_param_sigma_s(i);
    modelParams.a                   = gt_sigma_ori_scale*opt_param_sigma_s(i);
    modelParams.biasAmp             = opt_param_bias;
    modelParams.shape               = opt_param_shape;
    modelParams.scale               = opt_param_scale;
    modelParams.Cc                  = opt_param_Cc;
    modelParams.sigma_meta          = opt_param_sigma_meta;
    modelParams.guessRate           = opt_param_guessrate;
    
%     retData = getEstimatesPDFs(1:10:180, rvOriErr, modelParams);
    retData = getEstimationsPDF_cov(1:10:180, rvOriErr, modelParams);
    
    subplot(2, n_uncertainty_levels/2, i)
    hold on
    
    grpOriErr = resp_err_all_reshaped(i, :);
    histogram(grpOriErr, rvOriErr, Normalization="pdf");
    plot(rvOriErr, retData.analyticalPDF, LineWidth=1.5);
    
    xlabel("Orientation (deg)")
    ylabel("count")
    
    hold off
end

%% Plot results

figure

mean_sigma_s_stim = mean(sigma_s_stim, 2);

x = mean(resp_err_all_reshaped, 2);
x_m = median(resp_err_all_reshaped, 2);
y = std(resp_err_all_reshaped, 0, 2);
y_m = mad(resp_err_all_reshaped, 1, 2);

HC_idx = confidence_report_all_reshaped == 1;
LC_idx = confidence_report_all_reshaped == 0;

resp_HC = resp_err_all_reshaped;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all_reshaped;
resp_LC(~LC_idx) = NaN;

x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');
y_HC_m = mad(resp_HC, 1, 2);

x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');
y_LC_m = mad(resp_LC, 1, 2);

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
xlabel("\sigma_s(s) (sensory noise)")
ylabel("\sigma_m(s) (measurement noise)")
hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(mean_sigma_s_stim, y_HC, "filled", DisplayName="High confidence");
hold on
plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_HC, LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s_stim, y_LC, "filled", DisplayName="Low confidence");
plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_LC, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (sensory noise)")
ylabel("\sigma_m(s) (measurement noise)")
legend
hold off

subplot(2, 3, 4)

% Behavioral variability
scatter(mean_sigma_s_stim, y_m, "filled");
hold on
plot(mean_sigma_s_stim, anlytcl_mad_m_stim, LineWidth=1.5);
xlabel("\sigma_s (sensory noise)")
ylabel("MAD (measurement)")
hold off

subplot(2, 3, 5)

% Behavioral variability
scatter(mean_sigma_s_stim, y_HC_m, "filled", DisplayName="High confidence");
hold on
plot(mean_sigma_s_stim, anlytcl_mad_m_stim_HC, LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s_stim, y_LC_m, "filled", DisplayName="Low confidence");
plot(mean_sigma_s_stim, anlytcl_mad_m_stim_LC, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (sensory noise)")
ylabel("MAD (measurement)")
legend
hold off
