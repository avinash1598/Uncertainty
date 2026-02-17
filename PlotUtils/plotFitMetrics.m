function plotFitMetrics(cr_data, data)

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')

grpOriErr            = data.resp_err_all; 
confReport           = data.confidence_report_all;
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

trainFitAIC = zeros(K*nPerm, 1);
trainFitBIC = zeros(K*nPerm, 1);
testDataAIC = zeros(K*nPerm, 1);
testDataBIC = zeros(K*nPerm, 1);

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
        
        nll_ind_bestfit          = val;
        minfvalsCov(K*(h-1) + k) = val;
        params                   = fitParamsCovModel(idx, :);
        nllCov                   = minimizeError_cov(params, metaData);
        nllCovModel(K*(h-1) + k) = nllCov;
        
        fvalsCov( ( ( K*(h-1) + k - 1)*itr + 1) : (K*(h-1) + k)*itr ) = fvalsCovModel;
        
        % ind model
        fvalsIndModel     = Data_.resultsListInd{K*(h-1) + k}.f;
        fitParamsIndModel = Data_.resultsListInd{K*(h-1) + k}.x;
        [val, idx]          = min(fvalsIndModel);
        
        nll_cov_bestfit          = val;
        minfvalsInd(K*(h-1) + k) = val;
        params                   = fitParamsIndModel(idx, :);
        nllInd                   = minimizeError(params, metaData);
        nllIndModel(K*(h-1) + k) = nllInd;
        
        fvalsInd( ( (K*(h-1) + k - 1)*itr + 1) : (K*(h-1) + k)*itr ) = fvalsIndModel;
        
        deltaNLL(K*(h-1) + k) = nllInd - nllCov; % negative value better for cov data - this is better
        
        % train dataset AIC and BIC
        AIC = 2*10 + 2*nll_ind_bestfit - ( 2*9 + 2*nll_cov_bestfit); % abs val greater than 10 good
        BIC = 10*log(numel(trlErrors)*(K-1)/K) + 2*nll_ind_bestfit - ( 9*log(numel(trlErrors)*(K-1)/K) + 2*nll_cov_bestfit );
        
        trainFitAIC( K*(h-1) + k ) = AIC;
        trainFitBIC( K*(h-1) + k ) = BIC;

        % test dataset AIC and BIC
        AIC = 2*10 + 2*nllInd - ( 2*9 + 2*nllCov); % abs val greater than 10 good
        BIC = 10*log(numel(trlErrors)/K) + 2*nllInd - ( 9*log(numel(trlErrors)/K) + 2*nllCov );
        
        testDataAIC( K*(h-1) + k ) = AIC;
        testDataBIC( K*(h-1) + k ) = BIC;
    end
end

[~, pval] = ttest(testDataAIC)
[~, pval] = ttest(testDataBIC)

% All metrics on test dataset
figure
subplot(2, 2, 1)
hold on
bins = 1.7:0.01:1.9;
histogram(nllCovModel./1000, bins, DisplayName='cov model')
histogram(nllIndModel./1000, bins, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("Fit on test data")

subplot(2, 2, 2)
histogram(deltaNLL) 
ylabel("Count")
xlabel("delta NLL (ind - cov)")
legend
title("Fit on test data")

% [~, pval] = ttest(deltaNLL)

% This is probably not appropriate
% subplot(2, 2, 3)
% histogram(testDataAIC)
% ylabel("Count")
% xlabel("AIC")
% title("AIC (test dataset)")
% 
% subplot(2, 2, 4)
% histogram(testDataBIC)
% ylabel("Count")
% xlabel("BIC")
% title("BIC (test dataset)")

% All metrics on train dataset
figure

subplot(2, 2, 1)
hold on
histogram(fvalsCov./1000, DisplayName='cov model')
histogram(fvalsInd./1000, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("fvals (train dataset)")

subplot(2, 2, 2)
hold on
bins = 7:0.01:7.2;
histogram(minfvalsCov./1000, bins, DisplayName='cov model')
histogram(minfvalsInd./1000, bins, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("Bestfit fvals (train dataset)")

subplot(2, 2, 3)
histogram(trainFitAIC)
ylabel("Count")
xlabel("AIC")
title("AIC (train dataset)")

subplot(2, 2, 4)
histogram(trainFitBIC)
ylabel("Count")
xlabel("BIC")
title("BIC (train dataset)")

end


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

