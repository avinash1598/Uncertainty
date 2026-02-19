function cv_result = NLLCrossValidate(data, errBins, K, nPerm)

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/OptimizationUtils/')

if nargin < 3
    K = 5;
    nPerm = 6;
end

% Cross-validation
trlData   = convertToTrialData(data);
grpOriErr = trlData.grpOriErr;
N         = numel(grpOriErr(:));

% For each fold save nll test and train data
resultsListCov = cell(nPerm, K);
resultsListInd = cell(nPerm, K);
foldIDs        = cell(nPerm,1);
perms          = cell(nPerm,1);

parpool;

% Run k-fold cross-validation
parfor h=1:nPerm
    perm = randperm(N);
    foldID = mod(perm-1, K) + 1;
    
    foldIDs{h} = foldID;
    perms{h}   = perm;

    for k = 1:K
        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        fprintf( 'cross validation itr: %d \n', K*(h-1) + k) 
        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        
        % ---- Split trials ----
        testIdx  = (foldID == k);
        trainIdx = ~testIdx;
        
        % Run multi-start optimization for cov model
        results = Optimize(data, errBins, "cov", trainIdx);
        % resultsListCov{idx} = results;
        resultsListCov{h, k} = results;
        
        % Run multi-start optimization for ind model
        results = Optimize(data, errBins, "ind", trainIdx);
        % resultsListInd{idx} = results;
        resultsListInd{h, k} = results;
        
        % Should i pick the max p-val? or something else like mode? Look at
        % distribution maybe and then decide.
    end
end

Data.resultsListCov = resultsListCov;
Data.resultsListInd = resultsListInd;
Data.perms          = perms;
Data.foldIDs        = foldIDs;

cv_result = Data;

end
