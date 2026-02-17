close all
clear all

addpath('LL_scripts/')

data = load('COR15.mat'); % Changed reward function - affects perceptual variance - more relaxed
% data = load('../Data/COR16.mat'); % orientation dependent reward
% data = load('../Data/COR17.mat'); % ori dependent reward - changed reward function c1 = 5
% COR18 - c1=1, c2=0.3, -3 (all HC)
% COR19 - c1=1, c2=0.3, -0.5 (all HC)

stimOri = data.dat.stimOri;
reportedOri = data.dat.reportedOri;

rawError = reportedOri - stimOri;
rawOriError = mod(rawError + 90, 180) - 90;

median_val = median(rawOriError);
mad_val = median(abs(rawOriError - median_val));
mad_val = mad_val / 0.6745;
lower_bound = median_val - 3*mad_val;
upper_bound = median_val + 3*mad_val;
fltIdx = rawOriError < lower_bound | rawOriError > upper_bound;
rawOriError(fltIdx) = NaN;

data.dat.rawOriError = rawOriError;

selectedRawOriError = data.dat.rawOriError;

%% Prepare data structures
% TODO: Arrange the data in following format nlevels, norientations, ntrialPerOri
% Once I have this datastructure I can just use the existing analysis
% script

[G, contrastLevels, spreadLevels, stimDur] = findgroups(data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur);
grpIdxes             = unique(G);
n_uncertainty_levels = numel(grpIdxes);
nTrials              = numel(rawOriError(G == grpIdxes(1)));

% Arrange groups in increasing uncertainty level
fitStds = zeros(1, n_uncertainty_levels);

for i=1:n_uncertainty_levels
    grpIdx          = grpIdxes(i);
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    % grpRawErr_Flt   = data.dat.rawOriError(G == grpIdx); % rawOriErrorFlt
    
    fltOriErr = grpRawErr_Flt;
    pd = fitdist(fltOriErr(~isnan(fltOriErr)), 'Normal');
    
    mu = pd.mu;
    sigma = pd.sigma;
    fitStds(i) = sigma;
end

[B_, idx] = sort(fitStds);
sidx = grpIdxes(idx);

% Arrange in structured format - same as the one used for analysis
theta_true_all        = zeros(n_uncertainty_levels, nTrials);
theta_resp_all        = zeros(n_uncertainty_levels, nTrials); % Recorded theta based on user response
confidence_report_all = zeros(n_uncertainty_levels, nTrials);
resp_err_all          = zeros(n_uncertainty_levels, nTrials);

for i=1:n_uncertainty_levels
    grpIdx          = sidx(i);
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    grpStimOri      = data.dat.stimOri(G == grpIdx);
    grpOriResp      = data.dat.reportedOri(G == grpIdx);
    grpReportedConf = data.dat.reportedConf(G == grpIdx);
    
    theta_true_all(i, :)          = grpStimOri;
    theta_resp_all(i, :)          = grpOriResp;
    confidence_report_all(i, :)   = grpReportedConf;
    resp_err_all(i, :)            = grpRawErr_Flt;
end

grpOriErr = resp_err_all;

%% Get PDF for HC and LC
rvOriErr     = -90:3:90;

% Get PDFs from data for HC and LC
pdf_stim_LC = zeros( n_uncertainty_levels, numel(rvOriErr) );
pdf_stim_HC = zeros( n_uncertainty_levels, numel(rvOriErr) );

for i=1:n_uncertainty_levels

    cR = confidence_report_all(i, :);
    dataHC = resp_err_all(i, cR == 1);
    dataLC = resp_err_all(i, cR == 0);
    
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

HC_idx = confidence_report_all == 1;
LC_idx = confidence_report_all == 0;

resp_HC = resp_err_all;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all;
resp_LC(~LC_idx) = NaN;

std_HC = std(resp_HC, 0, 2, 'omitnan');
std_LC = std(resp_LC, 0, 2, 'omitnan');

metaData.rvOriErr     = rvOriErr;
metaData.pdf_stim_HC  = pdf_stim_HC;
metaData.pdf_stim_LC  = pdf_stim_LC;
metaData.targetStds   = std( resp_err_all, [], 2, "omitnan")';
metaData.std_HC       = std_HC';
metaData.std_LC       = std_LC';


%% Model comparison - fit both independent and cov model
nItrs = 100;

fValsCov = zeros(1, nItrs);
fValsInd = zeros(1, nItrs);

for i = 1:nItrs

disp(i)

success = false;
    
while ~success

try

    %%%%%%%%%%%%%%%%%%%%% Independent model
    param_sigma_s        = rand(1, n_uncertainty_levels); % Choose b such that average noise level ranges from low to high (relative to internal noise level)
    param_shape          = rand;
    param_scale          = rand;
    param_sigma_meta     = rand;
    param_Cc             = rand; 
    
    params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];

    % Fit model - independent model
    flagCov = false;
    
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
    
    if exitflag < 1
        error("Exit flag less than 1")
    end
    
    fValsInd(i) = fval;

    success = true;

catch ME
    fprintf('Error on iteration %d: %s\n', i, ME.message);
end

end

success = false;
    
while ~success

try

    %%%%%%%%%%%%%%%%%%%%% Covarying model
    param_sigma_s        = rand(1, n_uncertainty_levels); % Choose b such that average noise level ranges from low to high (relative to internal noise level)
    param_scale          = rand;
    param_sigma_meta     = rand;
    param_Cc             = rand; 
    
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
    
    if exitflag < 1
        error("Exit flag less than 1")
    end
    
    fValsCov(i) = fval;
    
    success = true;

catch ME
    fprintf('Error on iteration %d: %s\n', i, ME.message);
end

end

end

data.fValsCov = fValsCov;
data.fValsInd = fValsInd;
save("FltExpDataModelComparison.mat", "data")

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


