%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')

orientations     = 0:15:175; %linspace(0, 180, 18); %0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 250; %250;
b                = linspace(1, 2.2, 6); % linspace(1, 2.2, 8); Note: different minimum noise level (0.1). Choose b such that average noise level ranges from low to high (relative to internal noise level)
a                = 0.67.*b; %0.67.*b;   % Does a depend upon b? Yes
biasAmp          = 0.5; %0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% scale            = 0.5; %0.5;
% sigma_meta       = 0.6;
% Cc               = 0.5; 
% motor_noise_sigma = 0;
% 
% % biasAmp          = 0;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% % scale            = 343.3225;
% % sigma_meta       = 10.83;
% % Cc               = 0.05; 
% guessRate        = 0.1; %0.1; % While fitting try keeping it below 0.1 % For each trial with this prob sample uniformly from 0 to 179

% % a = [6.0984 8.2550 10.6963 8.8310 14.1120 16.5478 ];
% b = [11.3673   16.5831   17.3134   24.4308   30.7070   34.6668]; %[8.9120    8.3244   12.8387   10.8964   15.1622   18.2190 ]; %[8.7344 8.4603 12.9661 11.1144 14.5500 17.8371 ];
% % b = [6.9007 22.1322 21.0046 27.8742 36.0335 39.5992];
% a = 0.3186.*b; %0.8472
% biasAmp          = 1.6214; %9.9020 ; %0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
scale            = 163.3806; %54.7337 ; %0.5;
sigma_meta       = 9.7355; %521.5158 ;
Cc               = 0.0461 ; %0.1042 ; 
guessRate        = 0.3808; %0.1790 ;
motor_noise_sigma = 0;

% In actual data correct for bias

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = 6; %numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response

confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s = [12.0081   18.7643   19.3554   25.0146   32.1753   35.2551]';
sigma_s_stim = repmat(sigma_s, [1 n_theta]);
% sigma_s_stim = b' + a'*(abs(sind(2*orientations)));
% bias         = biasAmp*sind(2*orientations); 

% sigma_s_stim = b' + a'*(abs(sind(orientations - 90)));
bias         = biasAmp*sind(2*orientations); % This remains the same
bias(:) = 0;

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
        
        % Multiplicative
        shape = sigma_s_stim(l, i).^2 / scale; % divide by scale so that mean is sigma_s
        gain = gamrnd(shape, scale, [trials 1]);
        sigma_m_stim = sqrt( gain );
        mean_m_stim = theta_true + bias(i);
        
        % TODO: take into account bias?
        % Wrap the angle at the plotting stage. Note: warapping should be
        % performed according to the true angle.
        theta_est = mean_m_stim + sigma_m_stim .* randn(trials, 1);
        theta_est = mod(theta_est, 180); % Since this is orientation, wrap the angle between 0 and 180
        
        % Simulat guess rate
        guess_tl_idx = randi([1 trials], floor( trials*guessRate ), 1);
        guessOris = 180*rand(numel(guess_tl_idx), 1);
        theta_est(guess_tl_idx) = guessOris;
        
        % Add motor noise to estimate of theta
        theta_est = theta_est + motor_noise_sigma.*randn(trials, 1);
        
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
% data.params.b                     = b;
% data.params.a                     = a;
% data.params.biasAmp               = biasAmp;
data.params.scale                 = scale;
data.params.sigma_meta            = sigma_meta;
data.params.Cc                    = Cc;
data.params.guessRate             = guessRate;

% save('modelContOriData_cov.mat', "data")

%% Get analytical solution
anlytcl_sigma_m_stim = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_LC = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_mad_m_stim_LC = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    rvOriErr = -90:2:90; %
    
    modelParams.b                   = b(i);
    modelParams.a                   = a(i);
    modelParams.sigma_s             = data.params.sigma_s_reduced_model(i); % It would fit nicely if it were opt
    modelParams.biasAmp             = biasAmp;
    modelParams.scale               = scale;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    modelParams.guessRate           = guessRate;
    modelParams.sigma_m_shape1      = 1;
    modelParams.sigma_m_shape2      = 90;
    modelParams.biasShape           = 2;
    

%     tic
%     [~] = getPDFs_cov(orientations, modelParams, true);
%     elapsed_time = toc;
%     disp(['Execution time: ', num2str(elapsed_time), ' seconds']);
    
    retData = getPDFs_cov_reduced(modelParams, false);

%     tic
%     retData = getEstimationsPDF_cov_reduced(rvOriErr, modelParams);
%     elapsed_time = toc;
%     disp(['Execution time::::: ', num2str(elapsed_time), ' seconds']);
    
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
    
    modelParams.b                   = b(i);
    modelParams.a                   = a(i);
    modelParams.sigma_s             = data.params.sigma_s_reduced_model(i);
    modelParams.biasAmp             = biasAmp;
    modelParams.scale               = scale;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    modelParams.guessRate           = guessRate;
    modelParams.sigma_m_shape1      = 1;
    modelParams.sigma_m_shape2      = 90;
    modelParams.biasShape           = 2;
    

    subplot(2, n_uncertainty_levels/2, i)
    hold on
    grpOriErr = resp_err_all_reshaped(i, :);
    histogram(grpOriErr, rvOriErr, Normalization="pdf"); % Normalization="pdf"
    retData = getPDFs_cov_reduced(modelParams, false);
    %     retData = getPDFs_cov(modelParams);

    plot(retData.rvOriErrs, retData.analyticalPDF, LineWidth=1.5);
    
%     for j=[0.1 1 3]
%         errbn = -90:j:90;
%         
%         modelParams.b                   = b(i);
%         modelParams.a                   = a(i);
%         modelParams.sigma_s             = data.params.sigma_s_reduced_model(i);
%         modelParams.biasAmp             = biasAmp;
%         modelParams.scale               = scale;
%         modelParams.Cc                  = Cc;
%         modelParams.sigma_meta          = sigma_meta;
%         modelParams.guessRate           = guessRate;
%         
%         retData = getEstimationsPDF_cov(orientations, errbn, modelParams, true);
%     %     retData = getEstimationsPDF_cov_reduced(rvOriErr, modelParams);
% 
%         plot(retData.rvOriErrs, retData.analyticalPDF, LineWidth=1.5, DisplayName=""+j);
%     end
    
    %xlim([-90, 90])
    xlabel("Orientation (deg)")
    ylabel("count")
    % legend
    
    hold off
end

% PDFs by confidence and uncertainty
% Histogram by confidence
% Plot PDFs
figure 

for i=1:n_uncertainty_levels

    modelParams.b                   = b(i);
    modelParams.a                   = a(i);
    modelParams.sigma_s             = data.params.sigma_s_reduced_model(i);
    modelParams.biasAmp             = biasAmp;
    modelParams.scale               = scale;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    modelParams.guessRate           = guessRate;
    modelParams.sigma_m_shape1      = 1;
    modelParams.sigma_m_shape2      = 90;
    modelParams.biasShape           = 2;
    
%     retData = getPDFs_cov(orientations, modelParams, false);
    retData = getPDFs_cov_reduced(modelParams, false);
    

    cR = confidence_report_all_reshaped(i, :);
    dataHC = resp_err_all_reshaped(i, cR == 1);
    dataLC = resp_err_all_reshaped(i, cR == 0);
    
    subplot(2, n_uncertainty_levels, i);
    hold on
    
    histogram(dataLC, rvOriErr, Normalization="pdf");
    plot(rvOriErr, retData.analyticalPDF_LC, HandleVisibility="off", LineWidth=1.5)
    
    hold off

    xline(0, LineStyle="--")
    xlabel("Error (deg)")
    ylabel("count")
    legend
    % title("")

    subplot(2, n_uncertainty_levels, n_uncertainty_levels + i);
    hold on;

    histogram(dataHC, rvOriErr, Normalization="pdf");
    plot(rvOriErr, retData.analyticalPDF_HC, HandleVisibility="off", LineWidth=1.5)
    
    xline(0, LineStyle="--")
    xlabel("Error (deg)")
    ylabel("count")
    % title("")

    legend
    hold off

end

%% Plot results

figure

% sigma_s_stim = b' + a'*(abs(sind(1*orientations - 90)));
% sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 );
% 
% mean_sigma_s_stim = sigma_s_reduced_model; %mean(sigma_s_stim, 2);
% [~, idx] = sort(mean_sigma_s_stim);

idx = 1:6;
% mean_sigma_s_stim = mean(sigma_s_stim, 2);
mean_sigma_s_stim = data.params.sigma_s_reduced_model;

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
errorbar(mean_sigma_s_stim(idx), ...
    x(idx), y(idx), 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("\sigma_s(s)")
ylabel("Error")

subplot(2, 3, 2)

% Behavioral variability
scatter(mean_sigma_s_stim(idx), y, "filled");
hold on
plot(mean_sigma_s_stim(idx), anlytcl_sigma_m_stim(idx), LineWidth=1.5);
xlabel("\sigma_s(s) (sensory noise)")
ylabel("\sigma_m(s) (measurement noise)")
hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(mean_sigma_s_stim(idx), y_HC(idx), "filled", DisplayName="High confidence");
hold on
plot(mean_sigma_s_stim(idx), anlytcl_sigma_m_stim_HC(idx), LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s_stim(idx), y_LC(idx), "filled", DisplayName="Low confidence");
plot(mean_sigma_s_stim(idx), anlytcl_sigma_m_stim_LC(idx), LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (sensory noise)")
ylabel("\sigma_m(s) (measurement noise)")
legend
hold off

subplot(2, 3, 4)

% Behavioral variability
scatter(mean_sigma_s_stim(idx), y_m, "filled");
hold on
plot(mean_sigma_s_stim(idx), anlytcl_mad_m_stim(idx), LineWidth=1.5);
xlabel("\sigma_s (sensory noise)")
ylabel("MAD (measurement)")
hold off

subplot(2, 3, 5)

% Behavioral variability
scatter(mean_sigma_s_stim(idx), y_HC_m(idx), "filled", DisplayName="High confidence");
hold on
plot(mean_sigma_s_stim(idx), anlytcl_mad_m_stim_HC(idx), LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s_stim(idx), y_LC_m(idx), "filled", DisplayName="Low confidence");
plot(mean_sigma_s_stim(idx), anlytcl_mad_m_stim_LC(idx), LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (sensory noise)")
ylabel("MAD (measurement)")
legend
hold off