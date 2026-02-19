function nllData = computeNLL_CV(data, errBins, cv_data)

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')

% Exp data
trlData              = convertToTrialData(data);
n_uncertainty_levels = trlData.n_uncertainty_levels;
trlErrors            = trlData.trlErrors;
trlConfReports       = trlData.trlConfReports;
trlUncertaintyLevels = trlData.trlUncertaintyLevels;

warning("This is not computed on filtered data")
y_mad      = trlData.y_mad;
y_HC_mad   = trlData.y_HC_mad;
y_LC_mad   = trlData.y_LC_mad;

% CV params
nPerm   = numel( cv_data.foldIDs );
K       = numel( cv_data.resultsListCov ) / nPerm;
% n       = numel( cv_data.resultsListCov );
% itr     = numel( cv_data.resultsListCov{n}.f );

foldIDs = cv_data.foldIDs;

nllCovModel = []; %zeros(K*nPerm, 1);
nllIndModel = []; %zeros(K*nPerm, 1);
deltaNLL    = []; %zeros(K*nPerm, 1); % ind - cov

fvalsCov = []; % zeros(K*nPerm*30, 1);
fvalsInd = []; % zeros(K*nPerm*30, 1);

minfvalsCov = []; %zeros(K*nPerm, 1);
minfvalsInd = []; %zeros(K*nPerm, 1);

for h=1:nPerm
    foldID = foldIDs{h};
    
    for k = 1:K
        disp(K*(h-1) + k)

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
        
        % TODO: compute it on filtered data
        metaData.targetMADs    = y_mad;
        metaData.targetMADs_HC = y_HC_mad;
        metaData.targetMADS_LC = y_LC_mad;

        % Cov model
        fvalsCovModel     = cv_data.resultsListCov{h, k}.f;
        fitParamsCovModel = cv_data.resultsListCov{h, k}.x;
        [val, idx]        = min(fvalsCovModel);
        params            = fitParamsCovModel(idx, :);
        nllCov            = computeNLLCov(params, metaData);

        minfvalsCov = [minfvalsCov val];
        nllCovModel = [nllCovModel nllCov];
        fvalsCov    = [fvalsCov fvalsCovModel];
        % minfvalsCov(K*(h-1) + k) = val;
        % nllCovModel(K*(h-1) + k) = nllCov;
        % fvalsCov( ( ( K*(h-1) + k - 1)*itr + 1) : (K*(h-1) + k)*itr ) = fvalsCovModel;
        
        % ind model
        fvalsIndModel     = cv_data.resultsListInd{h, k}.f;
        fitParamsIndModel = cv_data.resultsListInd{h, k}.x;
        [val, idx]        = min(fvalsIndModel);
        params            = fitParamsIndModel(idx, :);
        nllInd            = computeNLL(params, metaData);
        
        minfvalsInd = [minfvalsInd val];
        nllIndModel = [nllIndModel nllInd];
        fvalsInd    = [fvalsInd fvalsIndModel];
        % minfvalsInd(K*(h-1) + k) = val;
        % nllIndModel(K*(h-1) + k) = nllInd;
        % fvalsInd( ( (K*(h-1) + k - 1)*itr + 1) : (K*(h-1) + k)*itr ) = fvalsIndModel;
        
        deltaNLL = [deltaNLL (nllCov - nllInd)];
        % deltaNLL(K*(h-1) + k) = nllCov - nllInd; % negative value better for cov data - this is better
    end
end

nllData.nllCovModel  = nllCovModel;
nllData.nllIndModel  = nllIndModel;
nllData.deltaNLL     = deltaNLL;

nllData.fvalsCov     = fvalsCov;
nllData.fvalsInd     = fvalsInd;

nllData.minfvalsCov  = minfvalsCov;
nllData.minfvalsInd  = minfvalsInd;

end
