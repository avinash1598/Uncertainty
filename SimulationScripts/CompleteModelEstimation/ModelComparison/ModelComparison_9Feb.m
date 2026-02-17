close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/Scripts/')

modelData = load('../modelContOriData.mat');

%% Get analytical solution for cov and independent model
cr_data = load("cross_validation_ind.mat");
nPerm   = numel( cr_data.Data_.foldIDs );
K       = numel( cr_data.Data_.resultsListCov ) / nPerm;
n       = numel( cr_data.Data_.resultsListCov );
itr     = numel( cr_data.Data_.resultsListCov{n}.f );

fValsCov  = zeros(n*itr,1);
fValsInd  = zeros(n*itr,1);
paramsCov = zeros(n*itr, 9);
paramsInd = zeros(n*itr, 10);

for i = 1:n
    fValsCov( (i - 1)*itr + 1: i*itr ) = cr_data.Data_.resultsListCov{i}.f;
    fValsInd( (i - 1)*itr + 1: i*itr ) = cr_data.Data_.resultsListInd{i}.f;

    paramsCov( (i - 1)*itr + 1: i*itr, : ) = cr_data.Data_.resultsListCov{i}.x;
    paramsInd( (i - 1)*itr + 1: i*itr, : ) = cr_data.Data_.resultsListInd{i}.x;
end

[~, idx] = min(fValsCov);
paramCovModel = paramsCov(idx, :);

[~, idx] = min(fValsInd);
paramIndModel = paramsInd(idx, :);

%% Plot fit metrics
plotFitMetrics(cr_data, modelData.data);

%% Plot result
% Cov
plotFitResult(modelData.data, paramCovModel, "cov");

% Ind
plotFitResult(modelData.data, paramIndModel, "ind");
