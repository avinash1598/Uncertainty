% TODO: multi start optimization but keep only the best fit - maybe do this
% during robustness check

close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')

modelData = load('../modelContOriData_cov.mat');
grpOriErr = modelData.data.err; 
confReport   = modelData.data.confReport;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr = reshape(grpOriErr, n_uncertainty_levels, []);
confReport = reshape(confReport, n_uncertainty_levels, []);
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

    [countHC, edges] = histcounts(dataHC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);

    [countLC, edges] = histcounts(dataLC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);
    
    binnedcount_reported_err_LC(i, :) = countLC;
    binnedcount_reported_err_HC(i, :) = countHC;
    
end

metaData.rvOriErr                     = rvOriErr;
metaData.binnedcount_reported_err_HC  = binnedcount_reported_err_HC;
metaData.binnedcount_reported_err_LC  = binnedcount_reported_err_LC;

AIC_values = zeros(1, 50);
BIC_values = zeros(1, 50);

for itr = 1:50

    success = false;

    while ~success
        try 
            % Fit cov model
            param_sigma_s        = std(grpOriErr, [], 2)';  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
            param_scale          = rand;
            param_sigma_meta     = rand;
            param_Cc             = rand; 
            
            params = [param_sigma_s param_scale param_sigma_meta param_Cc];
            nParams = numel(params); 
            
            % Objective function
            objFun = @(x) minimizeError_cov(x, grpOriErr, metaData);
            
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
            [optimalValues, fval, ~, ~] = fmincon(objFun, x0, ...
                [], [], [], [], lb, ub, [], options);
            
            disp('Optimal parameters:');
            disp(optimalValues);
            disp('Final objective value:');
            disp(fval / numel(grpOriErr(:))); % Per trial fval
            
            nll_cov = fval;

            success = true;

        catch ME
            disp(ME)
        end
    end
    
    success = false;

    while ~success
        try 

        % Fit independent model
        param_sigma_s        = std(grpOriErr, [], 2)';  % rand(1, n_uncertainty_levels); % Choose b such that average noise level ranges from low to high (relative to internal noise level)
        param_shape          = rand;
        param_scale          = rand;
        param_sigma_meta     = rand;
        param_Cc             = rand; 
        
        params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];
        nParams = numel(params); 
        
        % Objective function
        objFun = @(x) minimizeError(x, grpOriErr, metaData);
        
        % Bounds (ga requires finite bounds!)
        lb = zeros(size(params));     % same as before
        ub = []; % example finite upper bounds
        
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
        [optimalValues, fval, ~, ~] = fmincon(objFun, x0, ...
            [], [], [], [], lb, ub, [], options);
        
        disp('Optimal parameters:');
        disp(optimalValues);
        disp('Final objective value:');
        disp(fval / numel(grpOriErr(:))); % Per trial fval - nll
        
        nll_ind = fval;

        success = true;

        catch ME
            disp(ME)
        end
    end
    
    AIC = 2*10 + 2*nll_ind - ( 2*9 + 2*nll_cov);
    BIC = 10*log(numel(grpOriErr(:))) + 2*nll_ind - ( 9*log(numel(grpOriErr(:))) + 2*nll_cov );
    disp("Final AIC value (greater than 10 good)")
    disp(AIC)
    disp("Final BIC value (greater than 10 good)")
    disp(BIC)
    
    AIC_values(itr) = AIC;
    BIC_values(itr) = BIC;
end

ModelCmprsn_data.AIC = AIC_values;
ModelCmprsn_data.BIC = BIC_values;

save('ModelCmprsn_data_cov.mat', 'ModelCmprsn_data');

%% Plot result
load('ModelCmprsn_data_cov.mat')

figure
subplot(2, 2, 1)
histogram(ModelCmprsn_data.AIC, BinEdges=-100:10:100)
title("AIC")

subplot(2, 2, 2)
histogram(ModelCmprsn_data.BIC, BinEdges=-100:10:100)
title("BIC")


%% Loss function for optimization
function loss = minimizeError_cov(params, data, metaData)

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

%% Loss function for optimization
function loss = minimizeError(params, data, metaData)

nLevels = size(data, 1);

% Params
param_sigma_s        = params(1:nLevels);
param_shape          = params(nLevels + 1);
param_scale          = params(nLevels + 2);
param_sigma_meta     = params(nLevels + 3);
param_Cc             = params(nLevels + 4);

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
    modelParams.shape               = param_shape;   
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    
    retData = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    
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

