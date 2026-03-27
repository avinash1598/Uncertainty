close all
clear all
% ITR 191
warning("Make sure generative model is correct. While generating simulated " + ...
    "data guess rate should only be used for low confidence reports.")

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsTrialData\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

expData            = load('../Data/COR33.mat');     % Akash
% expData            = load('./Data/COR31.mat');   % Tien
% expData            = load('./Data/COR32.mat');   % Jiaming
% expData            = load('../Data/CORNFB02.mat');% Jonathan
% expData            = load('Data\CORNFB01.mat');  % Yichao
% For subject 6 => David: make sure to change the orientation dependent
% error

dataFileNames = ["COR31.mat" "COR33.mat" "CORNFB01.mat" "CORNFB02.mat" "CORNFB03.mat"];
fltSessionIdx = [0, 0, 0, 0, 2, 1];

uniqContrasts = unique( expData.dat.stimContrast );
uniqSpreads   = unique( expData.dat.stimSpread );
uniqDurations = unique( expData.dat.stimDur );

fltData = expData.dat( expData.dat.session > 0, :);
% fltData = fltData(fltData.stimSpread == uniqSpreads(1) & fltData.stimContrast == uniqContrasts(1), :);
fltData = fltData(fltData.stimDur == uniqDurations(2) & fltData.stimContrast == uniqContrasts(1), :);
f.dat = fltData; %data.dat; %fltData;
formattedData = formatExpData(f, false, true);

initCond             = getInitialConditions(formattedData);
n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now
%% Cov model

% ----------------------------------------
optParams.nStarts = 30; %300;
resultFullCov = Optimize(formattedData, initCond, "cov", [], optParams, "full");
% resultFullCov = Optimize(formattedData, initCond, "cov", [], optParams, "reduced");

% % Simulated data
% [~, idx] = min(resultFullCov.f);
% params   = resultFullCov.x(idx, :);
% % params(end - 2) = 0;
% data     = GenerateCovData(params, n_uncertainty_levels, 'full', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3));
% initCond = getInitialConditions(data);
% 
% optParams.nStarts = 100; %300
% resultFullCovSim = Optimize(data, initCond, "cov", [], optParams, 'full');

%% Ind model
% ----------------------------------------
optParams.nStarts = 30; %300;
resultFullInd = Optimize(formattedData, initCond, "ind", [], optParams, "full");
% resultFullInd = Optimize(formattedData, initCond, "ind", [], optParams, "reduced");

% % Simulated data
% [~, idx] = min(resultFullInd.f);
% params   = resultFullInd.x(idx, :);
% % params(end - 2) = 0;
% data     = GenerateIndData(params, n_uncertainty_levels, 'full', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3));
% initCond = getInitialConditions(data);
% 
% optParams.nStarts = 300; %300
% resultFullIndSim = Optimize(data, initCond, "ind", [], optParams, 'full');

%%
% load("J_Marginal_1.mat")

% resultFullCov = dataToSave.resultFullCov;
% resultFullInd = dataToSave.resultFullInd;

BinEdges = 21000:100:25000;

figure
hold on
histogram(resultFullCov.f, BinEdges, DisplayName="Cov"); %BinEdges,
histogram(resultFullInd.f, BinEdges, DisplayName="Ind");
%xlabel("fvals (log)")
xlabel("fvals")
ylabel("count")
legend
hold off

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

% Cov
errBins = -90:3:90;
[~, idx] = min(resultFullCov.f);
modelParamsCov = resultFullCov.x(idx, :);
plotFitResult(formattedData, modelParamsCov, initCond, "cov", "full", errBins)

% Ind
errBins = -90:3:90;
[~, idx] = min(resultFullInd.f);
modelParamsInd = resultFullInd.x(idx, :);
plotFitResult(formattedData, modelParamsInd, initCond, "ind", "full", errBins)


%%
% dataToSave.resultFullCov = resultFullCov;
% dataToSave.resultFullInd = resultFullInd;
% 
% save("J_Marginal_1.mat", "dataToSave");

