clear all
close all

addpath('LL_scripts/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

nTrialsList = [5 15 25 50 100 200];
nItr        = 30;

optRunData.fvals       = zeros(numel(nTrialsList), nItr);
optRunData.bias        = zeros(numel(nTrialsList), nItr);
optRunData.var         = zeros(numel(nTrialsList), nItr);
optRunData.gtParams    = zeros(numel(nTrialsList), nItr, 9);
optRunData.estParams   = zeros(numel(nTrialsList), nItr, 9);


for i = 1:numel(nTrialsList)
    disp("-----------------------------")
    nTrials = nTrialsList(i);
    
    % For this trial no get the 
    for j = 1:nItr
        disp(j)
        
        success = false;
        
        while ~success
        
            try

                [data, metaData]     = generateModelData(nTrials);
                grpOriErr            = metaData.grpOriErr;
                n_uncertainty_levels = metaData.n_uncertainty_levels;
                
                param_sigma_s        = std(grpOriErr, [], 2)';  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
                param_scale          = rand;
                param_sigma_meta     = rand;
                param_Cc             = rand; 
                
                params = [param_sigma_s param_scale param_sigma_meta param_Cc];
                
                % Fit model
                nParams = numel(params); 
                
                % Objective function
                objFun = @(x) minimizeError(x, grpOriErr, metaData);
                
                % Bounds (ga requires finite bounds!)
                lb = zeros(size(params));     % same as before
                ub = []; % example finite upper bounds
                
                % Fmincon
                % Optimization options for fmincon
                options = optimoptions('fmincon', ...
                    'Display', 'iter', ...
                    'Algorithm', 'sqp', ...          % or 'interior-point', 'trust-region-reflective', etc.
                    'MaxIterations', 1000, ...
                    'OptimalityTolerance', 1e-6, ...
                    'StepTolerance', 1e-6);
                
                % Initial guess (required for fmincon)
                x0 = params;   % start in the middle of bounds, for example
                
                % Run fmincon
                [optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
                    [], [], [], [], lb, ub, [], options);
                
                fval = fval / numel(grpOriErr(:));
                
                disp('Optimal parameters:');
                disp(optimalValues);
                disp('Final objective value:');
                disp(fval); % Per trial fval
                
                opt_param_sigma_s        = optimalValues(1:n_uncertainty_levels);
                opt_param_scale          = optimalValues(n_uncertainty_levels + 1);
                opt_param_sigma_meta     = optimalValues(n_uncertainty_levels + 2);
                opt_param_Cc             = optimalValues(n_uncertainty_levels + 3);
                
                % %% Goodness of parameter and curve fit
                % % Ground truth
                gt_sigma_s    = data.params.sigma_s_reduced_model;
                gt_scale      = data.params.scale;
                gt_sigma_meta = data.params.sigma_meta;
                gt_Cc         = data.params.Cc;
                
                % Display parameters
                for k =1:n_uncertainty_levels
                    fprintf("GT: %.4f, Fit: %.4f \n", gt_sigma_s(k), opt_param_sigma_s(k))
                end
                
                fprintf("GT: %.4f, Fit: %.4f \n", gt_scale, opt_param_scale)
                fprintf("GT: %.4f, Fit: %.4f \n", gt_sigma_meta, opt_param_sigma_meta)
                fprintf("GT: %.4f, Fit: %.4f \n", gt_Cc, opt_param_Cc)
        
                optRunData.fvals(i, j) = fval;
                optRunData.gtParams(i, j, :) = [gt_sigma_s gt_scale gt_sigma_meta gt_Cc];
                optRunData.estParams(i, j, :) = [opt_param_sigma_s opt_param_scale opt_param_sigma_meta opt_param_Cc];
                
                success = true;
                
            catch ME
                % Handle the error (log it, display, etc.)
                fprintf('Error on iteration %d: %s\n', i, ME.message);
                
                % You can pause or delay before retrying
                pause(1);

            end
        end
    end
    
    save('NLL_cov_bias_var_test_data.mat', "optRunData");
end

save('NLL_cov_bias_var_test_data.mat', "optRunData");
disp("done")

%% Calculate bias and variance
close all

% Distribution of fvals
figure

for i = 1:numel(nTrialsList)
    vals = optRunData.fvals(i, :);

    subplot(2, 3, i)
    histogram(vals, BinEdges=1:0.1:5);
    xlabel("fvals")
    ylabel("count")
    title(sprintf("NTrials %d", nTrialsList(i)*6*12));
end

idx = abs( optRunData.estParams - optRunData.gtParams ) > 1e1;
eP = optRunData.estParams; eP(idx) = 0;
gP = optRunData.gtParams; %gP(idx) = 0;
bias = squeeze( mean( eP, 2 ) - mean( gP, 2 ) );
variance = squeeze( var( ( eP - mean( eP, 2 ) ).^2, [], 2 ) ); 

% TODO: plot for each parameter
figure

for p = 1:9
    subplot(3, 3, p)
    hold on

    x = nTrialsList*6*12; %1:numel(nTrialsList);
    y1 = variance(:, p);
    y2 = bias(:, p).^2;

    yyaxis left
    plot(x, y1, 'LineWidth', 1.5, 'DisplayName', 'var')
    ylabel('Variance')
    
    yyaxis right
    plot(x, y2, 'LineWidth', 1.5, 'DisplayName', 'bias')
    ylabel('Bias')
    
    xlabel('Num trials')

%     plot(x, y1, DisplayName="var", LineWidth=1.5)
%     plot(x, y2, DisplayName="bias", LineWidth=1.5)

    hold off
%     xlabel("Num trials")
%     ylabel("Val")
    title(sprintf( "Parameter %d", p) );
end


%%
function [data, metaData] = generateModelData(ntrials_per_ori)

orientations     = 0:15:175; 
% ntrials_per_ori  = ntrials_per_ori; 
b                = linspace(1, 2.2, 6); 
a                = 0.67.*b; 
biasAmp          = 0.5;       
scale            = 0.5;
sigma_meta       = 0.6;
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
        
        % Multiplicative
        shape = sigma_s_stim(l, i).^2 / scale; % divide by scale so that mean is sigma_s
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
        theta_true_all(l, i, :)         = theta_true;
        theta_resp_all(l, i, :)         = theta_est;
        confidence_report_all(l, i, :)  = Vc > Cc;
    end
end

% Plot the performance curves
resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error

% resp_err_all_reshaped = reshape(resp_err_all, uncertainty_levels, []);
% confidence_report_all_reshaped = reshape(confidence_report_all, uncertainty_levels, []);

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

grpOriErr            = data.err; 
confReport           = data.confReport;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr    = reshape(grpOriErr, n_uncertainty_levels, []);
confReport   = reshape(confReport, n_uncertainty_levels, []);
rvOriErr     = -90:0.1:90;

% Get PDFs from data for HC and LC
binnedcount_reported_err_LC = zeros( n_uncertainty_levels, numel(rvOriErr) );
binnedcount_reported_err_HC = zeros( n_uncertainty_levels, numel(rvOriErr) );

for i=1:n_uncertainty_levels
    
    cR = confReport(i, :);
    dataHC = grpOriErr(i, cR == 1);
    dataLC = grpOriErr(i, cR == 0);
    
    centers = rvOriErr;
    binWidth = mean(diff(centers));
    edges = [centers - binWidth/2, centers(end) + binWidth/2];

    [countHC, ~] = histcounts(dataHC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);

    [countLC, ~] = histcounts(dataLC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);
    
    binnedcount_reported_err_LC(i, :) = countLC;
    binnedcount_reported_err_HC(i, :) = countHC;
    
end

metaData.n_uncertainty_levels         = n_uncertainty_levels;
metaData.grpOriErr                    = grpOriErr;
metaData.confReport                   = confReport;
metaData.rvOriErr                     = rvOriErr;
metaData.binnedcount_reported_err_HC  = binnedcount_reported_err_HC;
metaData.binnedcount_reported_err_LC  = binnedcount_reported_err_LC;

end


%% Loss function for optimization
function loss = minimizeError(params, data, metaData)

nLevels = size(data, 1);

% Params
param_sigma_s        = params(1:nLevels);
param_scale          = params(nLevels + 1);
param_sigma_meta     = params(nLevels + 2);
param_Cc             = params(nLevels + 3);

% Metadata
rvOriErr                     = metaData.rvOriErr;
binnedcount_reported_err_HC  = metaData.binnedcount_reported_err_HC;
binnedcount_reported_err_LC  = metaData.binnedcount_reported_err_LC;

currPdfFit_HC = zeros(nLevels, numel(rvOriErr));
currPdfFit_LC = zeros(nLevels, numel(rvOriErr));
curr_pHC       = zeros(nLevels, 1);
curr_pLC       = zeros(nLevels, 1);

for i=1:nLevels
    
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    
    retData = getEstimationsPDF_cov_reduced(rvOriErr, modelParams);
    
    % Data for NLL
    currPdfFit_HC(i, :) = retData.analyticalPDF_HC;
    currPdfFit_LC(i, :) = retData.analyticalPDF_LC;
    curr_pHC(i)         = retData.pHC;
    curr_pLC(i)         = retData.pLC;
    
end

% NLL loss
ll_HC = binnedcount_reported_err_HC .* log( currPdfFit_HC.*curr_pHC + eps );
ll_LC = binnedcount_reported_err_LC .* log( currPdfFit_LC.*curr_pLC + eps );

nll = - ( sum(ll_HC(:)) + sum(ll_LC(:)) );

loss = nll;

end

