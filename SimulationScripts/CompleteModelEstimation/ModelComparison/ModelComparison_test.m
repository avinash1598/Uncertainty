
%% fit HC and LC error distribution for each uncertainty level
close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')

modelData = load('modelContOriData.mat');
grpOriErr = modelData.data.err; 
confReport   = modelData.data.confReport;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr = reshape(grpOriErr, n_uncertainty_levels, []);
confReport = reshape(confReport, n_uncertainty_levels, []);
rvOriErr     = -90:0.1:90;

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


%% Model comparison - fit both independent and cov model


%%%%%%%%%%%%%%%%%%%%% Independent model
param_sigma_s        = rand(1, n_uncertainty_levels); % std(grpOriErr, [], 2)'; 
% param_shape          = rand;
param_scale          = rand;
param_sigma_meta     = rand;
param_Cc             = rand; 

% params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];
params = [param_sigma_s param_scale param_sigma_meta param_Cc];

% Fit model - independent model
flagCov = true;

% Objective function
objFun = @(x) minimizeError(x, grpOriErr, metaData, flagCov);

% Bounds (ga requires finite bounds!)
lb = zeros(size(params));     % same as before
ub = []; % example finite upper bounds

% Optimization options for fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...          % or 'interior-point', 'trust-region-reflective', etc.
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
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

% %%%%%%%%%%%%%%%%%%%%% Covarying model
% param_sigma_s        = rand(1, n_uncertainty_levels); % Choose b such that average noise level ranges from low to high (relative to internal noise level)
% param_scale          = rand;
% param_sigma_meta     = rand;
% param_Cc             = rand; 
% 
% params = [param_sigma_s param_scale param_sigma_meta param_Cc];
% 
% % Fit model - independent model
% flagCov = true;
% 
% % Objective function
% objFun = @(x) minimizeError(x, grpOriErr, metaData, flagCov);
% 
% % Bounds (ga requires finite bounds!)
% lb = zeros(size(params));     % same as before
% ub = []; % example finite upper bounds
% 
% % Optimization options for fmincon
% options = optimoptions('fmincon', ...
%     'Display', 'iter', ...
%     'Algorithm', 'sqp', ...          % or 'interior-point', 'trust-region-reflective', etc.
%     'MaxIterations', 1000, ...
%     'OptimalityTolerance', 1e-6, ...
%     'StepTolerance', 1e-6);
% 
% % Initial guess (required for fmincon)
% x0 = params;   % start in the middle of bounds, for example
% 
% % Run fmincon
% [optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
%     [], [], [], [], lb, ub, [], options);



%% Loss function for optimization

function loss = minimizeError(params, data, metaData, flagCov)

if flagCov
    loss = lossFnCovModel(params, data, metaData);
else
    loss = lossFnIndModel(params, data, metaData);
end

end

% Independent model
function loss = lossFnIndModel(params, data, metaData)

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

% Use common objective function
loss = sum( lossHC +  lossLC , 'all') + ...
    sum( ( curr_sigma_m - targetStds ).^2 ) + ...
    sum( ( curr_sigma_m_HC - targetStds_HC ).^2 ) + ...
    sum( ( curr_sigma_m_LC - targetStds_LC ).^2 );
    
end

% Cov model
function loss = lossFnCovModel(params, data, metaData)

nLevels = size(data, 1);

% Params
param_sigma_s        = params(1:nLevels);
param_scale          = params(nLevels + 1);
param_sigma_meta     = params(nLevels + 2);
param_Cc             = params(nLevels + 3);

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
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    
    retData = getEstimationsPDF_cov_reduced(rvOriErr, modelParams);
    
    currFit_HC(i, :) = retData.analyticalPDF_HC;
    currFit_LC(i, :) = retData.analyticalPDF_LC;
    curr_sigma_m(i)  = retData.E_sigma_m;
    curr_sigma_m_HC(i) = retData.E_sigma_m_HC;
    curr_sigma_m_LC(i) = retData.E_sigma_m_LC;
end

lossHC = ( currFit_HC - targetPDF_HC ).^2;
lossLC = ( currFit_LC - targetPDF_LC ).^2;


loss = sum( lossHC +  lossLC , 'all') + ...
    sum( ( curr_sigma_m - targetStds ).^2 ) + ...
    sum( ( curr_sigma_m_HC - targetStds_HC ).^2 ) + ...
    sum( ( curr_sigma_m_LC - targetStds_LC ).^2 );

end


