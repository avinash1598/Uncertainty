function result = Optimize(data, errBins, modelType, fltTrlIdx, optParams, fitType)

    assert(modelType == "cov" || modelType == "ind")

    if ~isempty(fitType)
        assert(fitType == "reduced" || fitType == "full")
    else
        fitType = "reduced";
    end

    if nargin < 5 || isempty(optParams)
        nStarts          = 20;
        hyperParamC1     = 100;
        randomGuessModel = true;
    else
        nStarts          = optParams.nStarts;
        hyperParamC1     = optParams.hyperParamC1;
        randomGuessModel = optParams.randomGuessModel;
    end
    
    fprintf("Params - nStarts: %d, hyperparam c1: %d", nStarts, hyperParamC1);

    addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')
    addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
    
    trlData              = convertToTrialData(data);
    grpOriErr            = trlData.grpOriErr;
    n_uncertainty_levels = trlData.n_uncertainty_levels;
    trlErrors            = trlData.trlErrors;
    trlConfReports       = trlData.trlConfReports;
    trlUncertaintyLevels = trlData.trlUncertaintyLevels;

    warning("This is not computed on filtered data. But maybe this is the right approach since this is th ground truth.")
    y_mad      = trlData.y_mad;
    y_HC_mad   = trlData.y_HC_mad;
    y_LC_mad   = trlData.y_LC_mad;
    
    if nargin > 3 && ~isempty(fltTrlIdx)
        trlErrors            = trlErrors(fltTrlIdx);
        trlConfReports       = trlConfReports(fltTrlIdx);
        trlUncertaintyLevels = trlUncertaintyLevels(fltTrlIdx);
    end
    
    % Compute mad on filtered data directly
%     metrics = computeMetricsFromTrlData(trlErrors, ...
%         trlConfReports, trlUncertaintyLevels);

    binnedData = buildBinnedData( ...
        n_uncertainty_levels, ...
        errBins, ...
        trlErrors, ...
        trlConfReports, ...
        trlUncertaintyLevels);
    
    metaData.n_levels      = n_uncertainty_levels;
    metaData.errBins       = errBins;
    metaData.binned_err_HC = binnedData.binned_err_HC;
    metaData.binned_err_LC = binnedData.binned_err_LC;
    metaData.targetMADs    = y_mad;
    metaData.targetMADs_HC = y_HC_mad;
    metaData.targetMADS_LC = y_LC_mad;
%     metaData.targetMADs    = metrics.mads;
%     metaData.targetMADs_HC = metrics.mads_HC;
%     metaData.targetMADS_LC = metrics.mads_LC;
    metaData.hyperParamC1  = hyperParamC1;
    metaData.randomGuessModel  = randomGuessModel;
    
    % Run multi-start optimization for cov model
    if fitType == "full"
        result = multiStartFitFull(grpOriErr, n_uncertainty_levels, metaData, modelType, nStarts);
    else
        result = multiStartFit(grpOriErr, n_uncertainty_levels, metaData, modelType, nStarts);
    end
    
end

%% Fucntions
% function to perform optimization (do with multistart option)
function results = multiStartFit(grpOriErr, n_uncertainty_levels, metaData, model, nStarts)
% model: cov, ind

if model == "cov"
    nParams = n_uncertainty_levels + 4;
elseif model == "ind"
    nParams = n_uncertainty_levels + 5;
else
    nParams = nan;
end

x_all = zeros(nStarts,nParams);
f_all = zeros(nStarts,1);

for itr = 1:nStarts

    fprintf( 'optimization itr: %d \n', itr) 
    success = false;
    
    while ~success
        try 
            
            param_sigma_s        = std(grpOriErr, [], 2)';  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
            param_scale          = rand;
            param_sigma_meta     = rand;
            param_Cc             = rand;
            param_guessrate      = 0.1*rand;
                
            if model == "cov"
                % params = [param_sigma_s param_scale param_sigma_meta param_Cc];
                params = [param_sigma_s param_scale param_sigma_meta param_Cc param_guessrate];
                objFun = @(x) computeNLLCov(x, metaData); % Objective function
            
            elseif model == "ind"
            
                param_shape = rand;
                % params      = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];
                params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc param_guessrate];
                objFun = @(x) computeNLL(x, metaData); % Objective function
            end
            
            % Bounds (ga requires finite bounds!)
            lb = zeros(size(params));     % same as before
            ub = inf( 1, numel(params) ); % example finite upper bounds
            ub(end) = 0.1; % Upper bound for last parameter i.e. guessrate
            
            warning('off','all')
            
            options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...          
                'MaxIterations', 1000, ...
                'MaxFunctionEvaluations', 20000);
            
            x0 = params;   % Initial guess (required for fmincon)
            
            [optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
                [], [], [], [], lb, ub, [], options);
            
            warning('on','all')
            
            disp(exitflag)
            disp(output.firstorderopt)
            
            if exitflag <= 0
                disp("fminconn failed")
                error('fminconn failed: %s', output.message)
            end

            x_all(itr, :) = optimalValues;
            f_all(itr)    = fval;
    
            success = true;

        catch ME
            disp(ME)
        end
    end

end

results.x = x_all;
results.f = f_all;

% First verify and then later pick the minimum nll
end


function results = multiStartFitFull(grpOriErr, n_uncertainty_levels, metaData, model, nStarts)
% model: cov, ind

% nStarts = 20;

if model == "cov"
    nParams = n_uncertainty_levels + 4 + 2;
elseif model == "ind"
    nParams = n_uncertainty_levels + 5 + 2;
else
    nParams = nan;
end

x_all = zeros(nStarts,nParams);
f_all = zeros(nStarts,1);

for itr = 1:nStarts

    fprintf( 'optimization itr: %d \n', itr) 
    success = false;
    
    while ~success
        try 
            
            param_sigma_s         = std(grpOriErr, [], 2)';  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
            param_scale           = rand;
            param_sigma_meta      = rand;
            param_Cc              = rand;
            param_guessrate       = 0.1*rand;
            param_sigma_ori_scale = rand;
            param_bias            = rand;
                
            if model == "cov"
                % params = [param_sigma_s param_scale param_sigma_meta param_Cc];
                params = [param_sigma_s param_scale param_sigma_meta param_Cc param_guessrate param_sigma_ori_scale param_bias];
                objFun = @(x) computeNLLCov(x, metaData, 'full'); % Objective function
            
            elseif model == "ind"
            
                param_shape = rand;
                % params      = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];
                params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc param_guessrate param_sigma_ori_scale param_bias];
                objFun = @(x) computeNLL(x, metaData, 'full'); % Objective function
            end
            
            % Bounds (ga requires finite bounds!)
            lb = zeros(size(params));     % same as before
            ub = inf( 1, numel(params) ); % example finite upper bounds
            ub(end - 2) = 0.1; % Upper bound for last parameter i.e. guessrate
            
            warning('off','all')
            
            options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...          
                'MaxIterations', 1000, ...
                'MaxFunctionEvaluations', 20000);
            
            x0 = params;   % Initial guess (required for fmincon)
            
            [optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
                [], [], [], [], lb, ub, [], options);
            
            warning('on','all')
            
            disp(exitflag)
            disp(output.firstorderopt)
            
            if exitflag <= 0
                disp("fminconn failed")
                error('fminconn failed: %s', output.message)
            end

            x_all(itr, :) = optimalValues;
            f_all(itr)    = fval;
    
            success = true;

        catch ME
            disp(ME)
        end
    end

end

results.x = x_all;
results.f = f_all;

% First verify and then later pick the minimum nll
end


