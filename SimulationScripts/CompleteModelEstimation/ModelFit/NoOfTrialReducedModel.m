clear all
close all

addpath('LL_scripts/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

% ntrials_per_oris=[25, 50, 100, 200, 500, 1000];
nItr = 50;
params_fitGoodnessData = zeros(1, nItr);
curve_fitGoodnessData  = zeros(1, nItr);

T = table([], [], [], [], ...
    [], [], [], [], ...
    [], [], [], [], ...
    [], [], [], [], ...
    'VariableNames', {'NumTrials', 'itr', 'sigma_s_1', 'sigma_s_2', 'sigma_s_3', ...
    'sigma_s_4', 'sigma_s_5', 'sigma_s_6', 'sigma_s_7', 'sigma_s_8', ...
    'shape', 'scale', 'sigma_meta', 'Cc', 'paramsFitLoss', 'curveFitLoss'});

dataToSave = {};

% for ntpo = numel(ntrials_per_oris)

for itr = 1:nItr
    
    disp(itr)
    
    success = false;
    
    while ~success

        try

            orientations     = linspace(0, 179, 10); %0:10:180; % 
            ntrials_per_ori  = 25;
            b                = linspace(0.01, 1.2, 8); % Choose b such that average noise level ranges from low to high (relative to internal noise level)
            a                = 0.67.*b;   % Does a depend upon b? Yes
            biasAmp          = 0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
            shape            = 2;
            scale            = 0.5;
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
            
            % Save model data
            data.stimOri     = theta_true_all;
            data.reportedOri = theta_resp_all;
            data.err         = resp_err_all;
            data.confReport  = confidence_report_all;
            
            data.params.sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
            data.params.b                     = b;
            data.params.a                     = a;
            data.params.biasAmp               = biasAmp;
            data.params.shape                 = shape;
            data.params.scale                 = scale;
            data.params.sigma_meta            = sigma_meta;
            data.params.Cc                    = Cc;
            
            dataToSave.groundTruth = [data.params.sigma_s_reduced_model ...
                data.params.shape ...
                data.params.scale ...
                data.params.sigma_meta ...
                data.params.Cc];
            
            save('modelContOriData_trialCntTest.mat', "data")
            
        
            %% Start fitting
            modelData = load('modelContOriData_trialCntTest.mat');
            grpOriErr = modelData.data.err;
            confReport   = modelData.data.confReport;
            n_uncertainty_levels = size(grpOriErr, 1);
            
            grpOriErr = reshape(grpOriErr, n_uncertainty_levels, []);
            confReport = reshape(confReport, n_uncertainty_levels, []);
            rvOriErr     = -90:0.1:90;
            
            param_sigma_s        = std(grpOriErr, [], 2)';  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
            param_shape          = rand;
            param_scale          = rand;
            param_sigma_meta     = rand;
            param_Cc             = rand; 
            
            params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];
            
            % Get PDFs from data for HC and LC
            pdf_stim_LC = zeros( n_uncertainty_levels, numel(rvOriErr) );
            pdf_stim_HC = zeros( n_uncertainty_levels, numel(rvOriErr) );
            
            for i=1:n_uncertainty_levels
            
                cR = confReport(i, :);
                dataHC = grpOriErr(i, cR == 1);
                dataLC = grpOriErr(i, cR == 0);
                
                centers = rvOriErr;
                binWidth = mean(diff(centers));
                edges = [centers - binWidth/2, centers(end) + binWidth/2];
            
                [pdfHC, edges] = histcounts(dataHC, ...
                    'Normalization', 'pdf', ...
                    'BinEdges', edges);
            
                [pdfLC, edges] = histcounts(dataLC, ...
                    'Normalization', 'pdf', ...
                    'BinEdges', edges);
                
                pdf_stim_HC(i, :) = pdfHC;
                pdf_stim_LC(i, :) = pdfLC;
                
            end
            
            HC_idx = confReport == 1;
            LC_idx = confReport == 0;
            
            resp_HC = grpOriErr;
            resp_HC(~HC_idx) = NaN;
            
            resp_LC = grpOriErr;
            resp_LC(~LC_idx) = NaN;
            
            std_HC = std(resp_HC, 0, 2, 'omitnan');
            std_LC = std(resp_LC, 0, 2, 'omitnan');

            metaData.rvOriErr     = rvOriErr;
            metaData.pdf_stim_HC  = pdf_stim_HC;
            metaData.pdf_stim_LC  = pdf_stim_LC;
            metaData.targetStds   = std( grpOriErr, [], 2 )';
            metaData.std_HC       = std_HC';
            metaData.std_LC       = std_LC';
            
            % Fit model
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
            [optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
                [], [], [], [], lb, ub, [], options);
            
            disp('Optimal parameters:');
            disp(optimalValues);
            disp('Final objective value:');
            disp(fval);
            
            opt_param_sigma_s        = optimalValues(1:n_uncertainty_levels);
            opt_param_shape          = optimalValues(n_uncertainty_levels + 1);
            opt_param_scale          = optimalValues(n_uncertainty_levels + 2);
            opt_param_sigma_meta     = optimalValues(n_uncertainty_levels + 3);
            opt_param_Cc             = optimalValues(n_uncertainty_levels + 4);
        
            %% Goodness of parameter and curve fit
        
            % Ground truth
            gt_sigma_s = modelData.data.params.sigma_s_reduced_model;
            gt_shape = modelData.data.params.shape;
            gt_scale = modelData.data.params.scale;
            gt_sigma_meta = modelData.data.params.sigma_meta;
            gt_Cc = modelData.data.params.Cc;
            
            loss_param_fit = sum( (opt_param_sigma_s - gt_sigma_s).^2 ) + ...
                (opt_param_shape - gt_shape).^2 + ...
                (opt_param_scale - gt_scale).^2 + ...
                (opt_param_sigma_meta - gt_sigma_meta).^2 + ...
                (opt_param_Cc - gt_Cc).^2;
            
            loss_curve_fit = 0;
            
            % Curve fit
            for i=1:n_uncertainty_levels
            
                modelParams.sigma_s             = gt_sigma_s(i);
                modelParams.shape               = gt_shape;   
                modelParams.scale               = gt_scale;
                modelParams.Cc                  = gt_sigma_meta;
                modelParams.sigma_meta          = gt_Cc;
                
                retData_gt = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
                
                modelParams.sigma_s             = opt_param_sigma_s(i);
                modelParams.shape               = opt_param_shape;   
                modelParams.scale               = opt_param_scale;
                modelParams.Cc                  = opt_param_Cc;
                modelParams.sigma_meta          = opt_param_sigma_meta;
                
                retData_fit = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
                
                loss_curve_fit = loss_curve_fit + ( retData_fit.analyticalPDF_LC - retData_gt.analyticalPDF_LC ).^2;
            
                loss_curve_fit = loss_curve_fit + ( retData_fit.analyticalPDF_HC - retData_gt.analyticalPDF_HC ).^2;
            end
            
            loss_curve_fit = sqrt(mean(loss_curve_fit));
            loss_param_fit = sqrt(loss_param_fit);
    
            params_fitGoodnessData(itr) = loss_curve_fit;
            curve_fitGoodnessData(itr)  = loss_param_fit;
            
            % 'VariableNames', {'NumTrials', 'itr', 'sigma_s_1', 'sigma_s_2', 'sigma_s_3', ...
            % 'sigma_s_4', 'sigma_s_5', 'sigma_s_6', 'sigma_s_7', 'sigma_s_8', ...
            % 'shape', 'scale', 'sigma_meta', 'Cc', 'paramsFitLoss', 'curveFitLoss'});
            row = {ntrials_per_ori*n_theta, ...
                itr, ...
                opt_param_sigma_s(1),  ...
                opt_param_sigma_s(2),  ...
                opt_param_sigma_s(3),  ...
                opt_param_sigma_s(4),  ...
                opt_param_sigma_s(5),  ...
                opt_param_sigma_s(6),  ...
                opt_param_sigma_s(7),  ...
                opt_param_sigma_s(8),  ...
                opt_param_shape,  ...
                opt_param_scale,  ...
                opt_param_sigma_meta, ... 
                opt_param_Cc, ...
                loss_param_fit, ...
                loss_curve_fit};
    
            T = [T; row];
            
            dataToSave.fitInfo = T;
            
            save("ModelFitMetrics_common_objective.mat", "dataToSave");
            
            fprintf('Saving done.')

            success = true;
            
        catch ME
            % Handle the error (log it, display, etc.)
            fprintf('Error on iteration %d: %s\n', i, ME.message);
            
            % You can pause or delay before retrying
            pause(1);

        end
    end
end

% end

dataToSave.fitInfo = T;

save("ModelFitMetrics_common_objective.mat", "dataToSave");


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
rvOriErr             = metaData.rvOriErr;
targetPDF_HC         = metaData.pdf_stim_HC;
targetPDF_LC         = metaData.pdf_stim_LC;
targetStds           = metaData.targetStds;
targetStds_HC        = metaData.std_HC;
targetStds_LC        = metaData.std_LC;

currFit_HC = zeros(nLevels, numel(rvOriErr));
currFit_LC = zeros(nLevels, numel(rvOriErr));
curr_sigma_m = zeros(1, nLevels);
curr_sigma_m_HC = zeros(1, nLevels);
curr_sigma_m_LC = zeros(1, nLevels);

for i=1:nLevels
    
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.shape               = param_shape;   
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    
    retData = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    
    currFit_HC(i, :) = retData.analyticalPDF_HC;
    currFit_LC(i, :) = retData.analyticalPDF_LC;
    curr_sigma_m(i)  = retData.E_sigma_m;
    curr_sigma_m_HC(i) = retData.E_sigma_m_HC;
    curr_sigma_m_LC(i) = retData.E_sigma_m_LC;
end

lossHC = ( currFit_HC - targetPDF_HC ).^2;
lossLC = ( currFit_LC - targetPDF_LC ).^2;

% loss = sum( lossHC + lossLC , 'all');

loss = sum( lossHC +  lossLC , 'all') + ...
    sum( ( curr_sigma_m - targetStds ).^2 ) + ...
    sum( ( curr_sigma_m_HC - targetStds_HC ).^2 ) + ...
    sum( ( curr_sigma_m_LC - targetStds_LC ).^2 );

% disp(loss)

end

