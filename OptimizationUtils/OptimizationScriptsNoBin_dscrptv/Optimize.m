function result = Optimize(data, initCond, modelType, fltTrlIdx, optParams, fitType)
    
    assert(modelType == "doublyStochastic" || modelType == "singlyStochastic")
    
    if ~isempty(fitType)
        assert(fitType == "reduced" || fitType == "full" || fitType == "full_jumbo")
    else
        fitType = "reduced";
    end
    
    if nargin < 5 || isempty(optParams)
        nStarts          = 20;
        %hyperParamC1     = 10; % Hyperparameter for data metrics
        %hyperParamC2     = 10; % Hyperparameter for oblique effect
        %randomGuessModel = true;
    else
        nStarts          = optParams.nStarts;
        %hyperParamC1     = optParams.hyperParamC1;
        %hyperParamC2     = optParams.hyperParamC2;
        %randomGuessModel = optParams.randomGuessModel;
    end
    
    fprintf('Params - nStarts: %d, Fittype: %s \n', nStarts, fitType);
    
    trlData              = convertToTrialData(data);
    grpOriErr            = trlData.grpOriErr;
    n_uncertainty_levels = trlData.n_uncertainty_levels;
    trlErrors            = trlData.trlErrors;
    trlConfReports       = trlData.trlConfReports;
    trlUncertaintyLevels = trlData.trlUncertaintyLevels;
    trlStimOris          = trlData.trlStimOris;
     
    if nargin > 3 && ~isempty(fltTrlIdx)
        trlErrors            = trlErrors(fltTrlIdx);
        trlConfReports       = trlConfReports(fltTrlIdx);
        trlUncertaintyLevels = trlUncertaintyLevels(fltTrlIdx);
        trlStimOris          = trlStimOris(fltTrlIdx);
    end
        
    metaData.n_levels             = n_uncertainty_levels;
    metaData.orientations         = data.orientations';
    metaData.trlErrors            = trlErrors;
    metaData.trlConfReports       = trlConfReports;
    metaData.trlUncertaintyLevels = trlUncertaintyLevels;
    metaData.trlStimOris          = trlStimOris;
    metaData.sigma_m_shape1       = initCond.sigma_m_shape1; % 1
    metaData.sigma_m_shape2       = initCond.sigma_m_shape2; % 90
    metaData.biasShape            = initCond.biasShape;      % 2
    metaData.a                    = initCond.a;
    metaData.biasAmp              = initCond.biasAmp;
    
    % Run multi-start optimization for cov model
    if fitType == "full"
        result = multiStartFitFull(n_uncertainty_levels, metaData, modelType, nStarts, initCond);
    else
        %result = multiStartFit(grpOriErr, n_uncertainty_levels, metaData, modelType, nStarts);
    end
    
end

%% Fucntions

function results = multiStartFitFull(n_uncertainty_levels, metaData, model, nStarts, initCond)
% model: cov, ind

% nStarts = 20;

if model == "doublyStochastic"
    nParams = 7;
elseif model == "singlyStochastic"
    nParams = 6;
else
    nParams = nan;
end

% Start pool only if none exists
p = gcp('nocreate');   
if isempty(p)
    parpool;           
end

x_all = zeros(nStarts,nParams);
f_all = zeros(nStarts,1);

% Initial conditions
b       = initCond.b;
a       = initCond.a;
biasAmp = initCond.biasAmp;

parfor itr = 1:nStarts
% for itr = 1:nStarts

    fprintf( 'optimization itr: %d \n', itr) 
    success = false;
    
    while ~success
        param_sigma_s         = rand; %rand(1, n_uncertainty_levels); %b;  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
        param_scale           = rand;
        param_sigma_meta      = rand;
        param_Cc              = rand;
        param_guessrate       = 0.01*rand;
        param_ori_scale       = rand; %a;
        param_bias            = biasAmp;
        
        if model == "doublyStochastic"
            params = [param_sigma_s ...
                param_scale ...
                param_sigma_meta ...
                param_Cc ...
                param_guessrate ...
                param_ori_scale ...
                param_bias];
            objFun = @(x) computeNLL_DoublyStochastic(x, metaData, 'full'); % Objective function
            
        elseif model == "singlyStochastic"
            params = [param_sigma_s ...
                param_sigma_meta ...
                param_Cc ...
                param_guessrate ...
                param_ori_scale ...
                param_bias];
            objFun = @(x) computeNLL_SinglyStochastic(x, metaData, 'full'); % Objective function
        end
        
        % Bounds (ga requires finite bounds!)
        lb = zeros(size(params));     % same as before
        ub = inf( 1, numel(params) ); % example finite upper bounds
        ub(end - 2) = 0.01; % literally set it to zero to see what happens
        
        warning('off','all')
        
        options = optimoptions('fmincon', ...
            'Display', 'iter', ...
            'Algorithm', 'sqp', ...          
            'MaxIterations', 1000, ...
            'MaxFunctionEvaluations', 20000, ...
            'OutputFcn', @stopAfterOneHour);
        
        x0 = params;   % Initial guess (required for fmincon)
        
        try 
        
            [optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
                [], [], [], [], lb, ub, [], options);
            
            disp(exitflag)
            disp(output.firstorderopt)
            
            if exitflag <= 0
                disp("fminconn failed")
                error('fminconn failed: %s', output.message)
            end

            success = true;

        catch ME
            %disp(ME)
            %keyboard
        end
        
        warning('on','all')
        
        if success
        
            x_all(itr, :) = optimalValues;
            f_all(itr)    = fval;

            disp("Itr complete ------------------------")
        end
   
    end

end

results.x = x_all;
results.f = f_all;

% First verify and then later pick the minimum nll
end

function stop = stopAfterOneHour(~, ~, state)

persistent startTime
stop = false;

switch state
    case 'init'
        startTime = tic;

    case 'iter'
        if toc(startTime) > (3000/2)
            stop = true;
            fprintf('Stopped: exceeded 0.5 hour.\n');
        end
end
end

