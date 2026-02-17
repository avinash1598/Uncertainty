close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Scripts')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/LLScriptsUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/OptimizationUtils/')

% expData            = load('Data\CORNFB01.mat'); % Yichao
expData            = load('./Data/COR33.mat'); % Akash
% expData            = load('Data\COR31.mat'); % Tien
% expData            = load('./Data/COR32.mat'); % Jiaming

fltData = expData.dat( expData.dat.session > 0 , :); 
f.dat = fltData; %data.dat; %fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

grpOriErr            = formattedData.resp_err_all; 
confReport           = formattedData.confidence_report_all;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr    = reshape(grpOriErr, n_uncertainty_levels, []);
confReport   = reshape(confReport, n_uncertainty_levels, []);
errBins      = -90:1:90;

levels = 1:6; % these are just labels
uncertainty_levels   = repmat(levels', [1 size(grpOriErr, 2)]); % Note: these are just labels and the values may not necessarily reflect actual uncertainty level

% Faltten data structures
trlErrors            = grpOriErr(:);
trlConfReports       = confReport(:);
trlUncertaintyLevels = uncertainty_levels(:);

%%
% Cross-validation
K = 5; % set to 15. 20 might be too much. Instead fix to 10 and do random perm.
N = numel(grpOriErr);
nPerm = 6;

% For each fold save nll test and train data
resultsListCov = cell(nPerm*K,1);
resultsListInd = cell(nPerm*K,1);
foldIDs        = cell(nPerm,1);
perms          = cell(nPerm,1);

% Run k-fold cross-validation
for h=1:nPerm
    perm = randperm(N);
    foldID = mod(perm-1, K) + 1;
    
    foldIDs{h} = foldID;
    perms{h}   = perm;

    for k = 1:K
        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        disp(k)
        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        
        % ---- Split trials ----
        testIdx  = (foldID == k);
        trainIdx = ~testIdx;
        
        % build binned data for train dataset
        binnedData = buildBinnedData( ...
            n_uncertainty_levels, ...
            errBins, ...
            trlErrors(trainIdx), ...
            trlConfReports(trainIdx), ...
            trlUncertaintyLevels(trainIdx));
        
        metaData.n_levels      = n_uncertainty_levels;
        metaData.errBins       = errBins;
        metaData.binned_err_HC = binnedData.binned_err_HC;
        metaData.binned_err_LC = binnedData.binned_err_LC;
        
        % Run multi-start optimization for cov model
        results = multiStartFit(grpOriErr, metaData, "cov");
        resultsListCov{K*(h-1) + k} = results;
        
        % Run multi-start optimization for ind model
        results = multiStartFit(grpOriErr, metaData, "ind");
        resultsListInd{K*(h-1) + k} = results;
        
        % Should i pick the max p-val? or something else like mode? Look at
        % distribution maybe and then decide.
    end
end

Data_.resultsListCov = resultsListCov;
Data_.resultsListInd = resultsListInd;
Data_.perms = perms;
Data_.foldIDs = foldIDs;
% save('cross_validation_tien.mat', 'Data_');

%% NLL of test dataset

cr_data = load('cross_validation_akash.mat');
nPerm   = numel( cr_data.Data_.foldIDs );
K       = numel( cr_data.Data_.resultsListCov ) / nPerm;
n       = numel( cr_data.Data_.resultsListCov );
itr     = numel( cr_data.Data_.resultsListCov{n}.f );

foldIDs = cr_data.Data_.foldIDs;
Data_   = cr_data.Data_;

nllCovModel = zeros(K*nPerm, 1);
nllIndModel = zeros(K*nPerm, 1);
deltaNLL    = zeros(K*nPerm, 1);

fvalsCov = zeros(K*nPerm*30, 1);
fvalsInd = zeros(K*nPerm*30, 1);

minfvalsCov = zeros(K*nPerm, 1);
minfvalsInd = zeros(K*nPerm, 1);


for h=1:nPerm
    foldID = foldIDs{h};
    
    for k = 1:K
        % ---- Split trials ----
        testIdx  = (foldID == k);
        
        % build binned data for train dataset
        binnedData = buildBinnedData( ...
            n_uncertainty_levels, ...
            errBins, ...
            trlErrors(testIdx), ...
            trlConfReports(testIdx), ...
            trlUncertaintyLevels(testIdx));
        
        metaData.n_levels      = n_uncertainty_levels;
        metaData.errBins       = errBins;
        metaData.binned_err_LC = binnedData.binned_err_LC;
        metaData.binned_err_HC = binnedData.binned_err_HC;
        
        % Cov model
        fvalsCovModel     = Data_.resultsListCov{K*(h-1) + k}.f;
        fitParamsCovModel = Data_.resultsListCov{K*(h-1) + k}.x;
        [val, idx]          = min(fvalsCovModel);
       
        minfvalsCov(K*(h-1) + k) = val;
        params                   = fitParamsCovModel(idx, :);
        nllCov                   = minimizeError_cov(params, metaData);
        nllCovModel(K*(h-1) + k) = nllCov;
        
        fvalsCov( ( ( K*(h-1) + k - 1)*30 + 1) : (K*(h-1) + k)*30 ) = fvalsCovModel;
        
        % ind model
        fvalsIndModel     = Data_.resultsListInd{K*(h-1) + k}.f;
        fitParamsIndModel = Data_.resultsListInd{K*(h-1) + k}.x;
        [val, idx]          = min(fvalsIndModel);
        
        minfvalsInd(K*(h-1) + k) = val;
        params                   = fitParamsIndModel(idx, :);
        nllInd                   = minimizeError(params, metaData);
        nllIndModel(K*(h-1) + k) = nllInd;
        
        fvalsInd( ( (K*(h-1) + k - 1)*30 + 1) : (K*(h-1) + k)*30 ) = fvalsIndModel;
        
        deltaNLL(K*(h-1) + k) = nllCov - nllInd; % negative value better for cov data - this is better
    end
end

figure
subplot(2, 2, 1)
hold on
% bins = 480:5:550;
histogram(nllCovModel, DisplayName='cov model')
histogram(nllIndModel, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("Fit on test data (GT: cov model data)")

subplot(2, 2, 2)
histogram(deltaNLL) %  -3:0.5:2
ylabel("Count")
xlabel("delta NLL (cov - ind)")
legend
title("Fit on test data (GT: cov model data)")

subplot(2, 2, 3)
hold on
bins = 3.3:0.01:3.8; %4.4:0.01:4.8;
histogram(fvalsCov./1000, DisplayName='cov model')
histogram(fvalsInd./1000, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("fvals (train dataset)")

subplot(2, 2, 4)
hold on
bins = 3.3:0.01:3.8; %4.4:0.01:4.8;
histogram(minfvalsCov./1000, DisplayName='cov model')
histogram(minfvalsInd./1000, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("bestfit fvals (train dataset)")


%% Fucntions
function binnedData = buildBinnedData(n_levels, errBins, trlErrs, trlConfReports, trlUncertaintyLevels)

% Get PDFs from data for HC and LC
binned_err_LC = zeros( n_levels, numel(errBins) );
binned_err_HC = zeros( n_levels, numel(errBins) );

for i=1:n_levels
    
    fltIdx = trlUncertaintyLevels == i;

    cR = trlConfReports(fltIdx);
    fltErr = trlErrs(fltIdx);

    dataHC = fltErr(cR == 1);
    dataLC = fltErr(cR == 0);
    
    centers = errBins;
    binWidth = mean(diff(centers));
    edges = [centers - binWidth/2, centers(end) + binWidth/2];
    
    [countHC, ~] = histcounts(dataHC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);

    [countLC, ~] = histcounts(dataLC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);
    
    binned_err_LC(i, :) = countLC;
    binned_err_HC(i, :) = countHC;
    
end

binnedData.binned_err_LC = binned_err_LC;
binnedData.binned_err_HC = binned_err_HC;

end

% function to perform optimization (do with multistart option)
function results = multiStartFit(grpOriErr, metaData, model)
% model: cov, ind

nStarts = 30;

if model == "cov"
    nParams = 9;
elseif model == "ind"
    nParams = 10;
else
    nParams = nan;
end

x_all = zeros(nStarts,nParams);
f_all = zeros(nStarts,1);

for itr = 1:nStarts

    success = false;

    while ~success
        try 

            param_sigma_s        = std(grpOriErr, [], 2)';  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
            param_scale          = rand;
            param_sigma_meta     = rand;
            param_Cc             = rand; 
                
            if model == "cov"
                params = [param_sigma_s param_scale param_sigma_meta param_Cc];
                objFun = @(x) minimizeError_cov(x, metaData); % Objective function
            
            elseif model == "ind"
            
                param_shape = rand;
                params      = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];
                objFun = @(x) minimizeError(x, metaData); % Objective function
            end
            
            % Bounds (ga requires finite bounds!)
            lb = zeros(size(params));     % same as before
            ub = [];                      % example finite upper bounds
            
            options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...          
                'MaxIterations', 1000, ...
                'MaxFunctionEvaluations', 20000, ...
                'OptimalityTolerance', 1e-6, ...
                'StepTolerance', 1e-6);
            
            x0 = params;   % Initial guess (required for fmincon)
            
            warning('off','all')

            [optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
                [], [], [], [], lb, ub, [], options);
            
            disp(exitflag)
            disp(output.firstorderopt)

            if exitflag <= 0
                error('fminconn failed: %s', output.message)
            end
            
            warning('on','all')
            
            x_all(itr, :) = optimalValues;
            f_all(itr)    = fval;
            
            success = true;

        catch ME
            % disp(ME)
            %disp("err")
        end
    end

end

results.x = x_all;
results.f = f_all;

% First verify and then later pick the minimum nll
end

% Loss function for optimization
function loss = minimizeError_cov(params, metaData)

nLevels = metaData.n_levels;

% Params
param_sigma_s        = params(1:nLevels);
param_scale          = params(nLevels + 1);
param_sigma_meta     = params(nLevels + 2);
param_Cc             = params(nLevels + 3);

% Metadata
errBins        = metaData.errBins;
binned_err_HC  = metaData.binned_err_HC;
binned_err_LC  = metaData.binned_err_LC;

currPdfFit_HC = zeros(nLevels, numel(errBins));
currPdfFit_LC = zeros(nLevels, numel(errBins));
curr_pHC       = zeros(nLevels, 1);
curr_pLC       = zeros(nLevels, 1);

for i=1:nLevels
    
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = 0;
    
    retData = getEstimationsPDF_cov_reduced(errBins, modelParams);
    
    % Data for NLL
    currPdfFit_HC(i, :) = retData.analyticalPDF_HC;
    currPdfFit_LC(i, :) = retData.analyticalPDF_LC;
    curr_pHC(i)         = retData.pHC;
    curr_pLC(i)         = retData.pLC;
    
end

% NLL loss
ll_HC = binned_err_HC .* log( currPdfFit_HC.*curr_pHC + eps );
ll_LC = binned_err_LC .* log( currPdfFit_LC.*curr_pLC + eps );

nll = - ( sum(ll_HC(:)) + sum(ll_LC(:)) );

loss = nll;

end

% Loss function for optimization
function loss = minimizeError(params, metaData)

nLevels = metaData.n_levels;

% Params
param_sigma_s        = params(1:nLevels);
param_shape          = params(nLevels + 1);
param_scale          = params(nLevels + 2);
param_sigma_meta     = params(nLevels + 3);
param_Cc             = params(nLevels + 4);

% Metadata
errBins        = metaData.errBins;
binned_err_HC  = metaData.binned_err_HC;
binned_err_LC  = metaData.binned_err_LC;

currPdfFit_HC = zeros(nLevels, numel(errBins));
currPdfFit_LC = zeros(nLevels, numel(errBins));
curr_pHC       = zeros(nLevels, 1);
curr_pLC       = zeros(nLevels, 1);

for i=1:nLevels
    
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.shape               = param_shape;   
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = 0;
    
    retData = getEstimatesPDFs_reduced_model(errBins, modelParams);
    
    % Data for NLL
    currPdfFit_HC(i, :) = retData.analyticalPDF_HC;
    currPdfFit_LC(i, :) = retData.analyticalPDF_LC;
    curr_pHC(i)         = retData.pHC;
    curr_pLC(i)         = retData.pLC;
    
end

% NLL loss
ll_HC = binned_err_HC .* log( currPdfFit_HC.*curr_pHC + eps );
ll_LC = binned_err_LC .* log( currPdfFit_LC.*curr_pLC + eps );

nll = - ( sum(ll_HC(:)) + sum(ll_LC(:)) );

loss = nll; 

end

