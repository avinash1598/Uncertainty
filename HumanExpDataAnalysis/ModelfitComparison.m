close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Scripts/')

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR31.mat'); % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR32.mat'); % Jiaming
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR33.mat');   % Akash
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/ExpScript/CORNFB01.mat');   % Yichao

%% Format data
fltData       = data.dat( data.dat.session > 0 , :); 
f.dat         = fltData; 
formattedData = formatExpData(f, false, false);  % Keep sortByMAD to false (it is set to false in NLL script)

%% Get analytical solution for cov and independent model
% cr_data = load("cross_validation_tien.mat");
cr_data = load("cross_validation_akash.mat");
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
plotFitMetrics(cr_data, formattedData); 

%% Plot result
% Cov
plotFitResult(formattedData, paramCovModel, "cov");

% Ind
plotFitResult(formattedData, paramIndModel, "ind");
