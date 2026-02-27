%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')

% orientations     = linspace(0, 179, 30); %0:10:180; % 
orientations     = linspace(0, 179, 10); % Alert!!!! This has impact on quality of analytical solution
ntrials_per_ori  = 2000; %1000;
b                = linspace(0.1, 1.5, 6); % 1.2 % Choose b such that average noise level ranges from low to high (relative to internal noise level)
a                = 0.67.*b; %0.67   % Does a depend upon b? Yes
% biasAmp          = 0; %10       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% shape            = 0.848;
% scale            = 338.1997;
% sigma_meta       = 41.5023;
% Cc               = 0.109; 
% guessRate        = 0;

biasAmp          = 2; %0.5; %0.5; % problem at 2       % Does bias depend upon uncertainty level? No. This bias level seems okay.
shape            = 2;
scale            = 0.5; %0.5;
sigma_meta       = 0.2;
Cc               = 0.7; %0.7
guessRate        = 0.1; %0.08;

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(b);
% n_uncertainty_levels     = numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s_stim = b' + a'*(abs(sind(( orientations - 90 ) )));
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

data.params.sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
data.params.b                     = b;
data.params.a                     = a;
data.params.biasAmp               = biasAmp;
data.params.shape                 = shape;
data.params.scale                 = scale;
data.params.sigma_meta            = sigma_meta;
data.params.Cc                    = Cc;
data.params.guessRate             = guessRate;

% save('modelContOriData.mat', "data")


%% Get analytical solution
anlytcl_sigma_m_stim = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_LC = zeros(1, uncertainty_levels);
anlytcl_mad_m_byOri   = zeros(uncertainty_levels, numel(orientations));
anlytcl_mad_m_stim    = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim_LC = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    rvOriErr = -90:0.5:90;
    
    modelParams.b                   = b(i);
    modelParams.a                   = a(i);
    modelParams.sigma_s             = data.params.sigma_s_reduced_model(i); % It would fit nicely if it were opt
    modelParams.biasAmp             = biasAmp;
    modelParams.shape               = shape;
    modelParams.scale               = scale;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    modelParams.guessRate           = guessRate;
    
    tic
    [~] = getEstimatesPDFs(orientations, rvOriErr, modelParams, true);
    %[~] = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    elapsed_time = toc;
    disp(['Execution time: ', num2str(elapsed_time), ' seconds']);
    
    retData = getEstimatesPDFs(orientations, rvOriErr, modelParams, false);
    % retData = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);

    anlytcl_sigma_m_stim(i)    = retData.E_sigma_m;
    anlytcl_sigma_m_stim_HC(i) = retData.E_sigma_m_HC;
    anlytcl_sigma_m_stim_LC(i) = retData.E_sigma_m_LC;
    
    anlytcl_mad_m_stim(i)      = retData.mad_m;
    anlytcl_mad_m_stim_HC(i)   = retData.mad_m_HC;
    anlytcl_mad_m_stim_LC(i)   = retData.mad_m_LC;
    anlytcl_mad_m_byOri(i, :)  = retData.mad_m_by_ori';
 
end

%%
% PDFs by uncertainty
figure

n_uncertainty_levels = numel(b);

for i=1:n_uncertainty_levels
    
    modelParams.b                   = b(i);
    modelParams.a                   = a(i);
    modelParams.sigma_s             = data.params.sigma_s_reduced_model(i);
    modelParams.biasAmp             = biasAmp;
    modelParams.shape               = shape;
    modelParams.scale               = scale;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    modelParams.guessRate           = guessRate;
    
    % retData = getEstimatesPDFs(1:10:180, rvOriErr, modelParams);
    retData = getEstimatesPDFs(orientations, rvOriErr, modelParams, false);
%     retData = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    
    subplot(3, n_uncertainty_levels, i)
    hold on
    
    grpOriErr = resp_err_all_reshaped(i, :);
    histogram(grpOriErr, rvOriErr, Normalization="pdf");
    plot(rvOriErr, retData.analyticalPDF, LineWidth=1.5);
    
    xlabel("Orientation (deg)")
    ylabel("count")
    
    hold off

    % HC
    subplot(3, n_uncertainty_levels, 6+i)
    hold on
    
    confReport_ = confidence_report_all_reshaped(i, :);
    
    errHC = grpOriErr(confReport_ == 1);
    histogram(errHC, rvOriErr, Normalization="pdf");
    plot(rvOriErr, retData.analyticalPDF_HC, LineWidth=1.5);
    
    xlabel("Orientation (deg)")
    ylabel("count")
    title("HC")
    
    hold off

    % LC
    subplot(3, n_uncertainty_levels, 12+i)
    hold on
    
    errHC = grpOriErr(confReport_ == 0);
    histogram(errHC, rvOriErr, Normalization="pdf");
    plot(rvOriErr, retData.analyticalPDF_LC, LineWidth=1.5);
    
    xlabel("Orientation (deg)")
    ylabel("count")
    title("LC")
    
    hold off
end

%% MAD by ori
figure

for i=1:n_uncertainty_levels
    
    x = squeeze( resp_err_all(i, :, :) );
    madByOri = mad(x, 1, 2);

    subplot(2, n_uncertainty_levels/2, i)
    hold on
    
    scatter(orientations, madByOri)
    plot(orientations, anlytcl_mad_m_byOri(i, :), LineWidth=1.5);
    
    xlabel("Orientation (deg)")
    ylabel("MAD")
    
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
