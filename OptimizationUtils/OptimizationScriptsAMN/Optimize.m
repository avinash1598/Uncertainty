function result = Optimize(data, initCond, fltTrlIdx, optParams)
    
    if nargin < 4 || isempty(optParams)
        nStarts          = 20;
    else
        nStarts          = optParams.nStarts;
    end
    
    fprintf('Params - nStarts: %d \n', nStarts);
    
    trlData              = convertToTrialData(data);
    grpOriErr            = trlData.grpOriErr;
    n_uncertainty_levels = trlData.n_uncertainty_levels;
    trlErrors            = trlData.trlErrors;
    trlConfReports       = trlData.trlConfReports;
    trlUncertaintyLevels = trlData.trlUncertaintyLevels;
    trlStimOris          = trlData.trlStimOris;
     
    if nargin > 2 && ~isempty(fltTrlIdx)
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
    result = multiStartFit(grpOriErr, n_uncertainty_levels, metaData, nStarts);
    
end

%% Fucntions
% function to perform optimization (do with multistart option)
function results = multiStartFit(grpOriErr, n_uncertainty_levels, metaData, nStarts)
% model: additive, multiplicative

nParams = 2*n_uncertainty_levels + 4;

% Start pool only if none exists
p = gcp('nocreate');   
if isempty(p)
    parpool;           
end

x_all = zeros(nStarts,nParams);
f_all = zeros(nStarts,1);

% for itr = 1:nStarts
parfor itr = 1:nStarts

    fprintf( 'optimization itr: %d \n', itr) 
    success = false;
    
    while ~success
        try 
            
            param_sigma_ext       = std(grpOriErr, [], 2)';  
            param_ext_noise_scale = rand(1, n_uncertainty_levels);
            param_scale           = rand;
            param_sigma_meta      = rand;
            param_Cc              = rand;
            param_guessrate       = 0.1*rand;
            
            params = [param_sigma_ext ...
                param_ext_noise_scale...
                param_scale...
                param_sigma_meta...
                param_Cc...
                param_guessrate];
            objFun = @(x) computeNLL(x, metaData); % Objective function
            
            % Bounds (ga requires finite bounds!)
            lb      = zeros(size(params));     % same as before
            ub      = inf( 1, numel(params) ); % example finite upper bounds
            ub(end) = 1; % Upper bound for last parameter i.e. guessrate
            
            warning('off','all')
            
            options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ... % % 'Algorithm', 'sqp', ...    interior-point
                'MaxIterations', 1000, ...
                'MaxFunctionEvaluations', 20000, ...
                'OutputFcn', @stopAfterOneHour);
            
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
            %disp(ME)
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
