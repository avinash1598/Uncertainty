close all
clear all
restoredefaultpath

% get(groot,"default")
set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',1, ...
    'defaultAxesBox','off', ...
    'defaultAxesTickDir','out', ...
    'defaultLineLineWidth',1.5, ...
    'defaultAxesLabelFontSizeMultiplier',1.1, ...
    'defaultAxesTitleFontWeight','normal', ...
    'defaultLegendBox','off');

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')
% 
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsTrialData\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

expData            = load('./Data/COR33.mat');     % Akash
% expData            = load('./Data/COR31.mat');   % Tien
% expData            = load('./Data/COR32.mat');   % Jiaming
% expData            = load('./Data/CORNFB02.mat');% Jonathan
% expData            = load('Data\CORNFB01.mat');  % Yichao
% expData            = load('Data\CORNFB01.mat');  % Yichao
% For subject 6 => David: make sure to change the orientation dependent
% error

fltSessionIdx = [0, 0, 0, 0, 2, 1];
fltData       = expData.dat( expData.dat.session > 0 , :);  % TODO: change session number
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

initCond      = getInitialConditions(formattedData);

n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now

%%
% load('./Data/ModelFitQualityTest_Akash.mat');
load('ModelFitQualityTest_Yichao.mat');
plotFitQuality(formattedData, dataToSave)
% load('./Data/ModelFitQualityTest_Jonathan.mat');

% Fit: guess rate is weirdly high!

% %% Plot fit results (Reduced)
% errBins = -90:3:90;
% [~, idx] = min(dataToSave.resultReducedCov.f);
% modelParams = dataToSave.resultReducedCov.x(idx, :);
% plotFitResult(formattedData, modelParams, initCond, "cov", "reduced", errBins)
% % NLL for each trial
% nllsData = getNLLForEachTrial(formattedData, initCond, "cov", "reduced", modelParams);
% 
% % Simulation - GR too high (red flag)
% % There is something about guess rate that impacts recoverability
% [~, idx] = min(dataToSave.resultReducedCovSim.f);
% modelParamsSim = dataToSave.resultReducedCovSim.x(idx, :);
% data = GenerateCovData( modelParamsSim, n_uncertainty_levels, 'reduced', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3) );
% initCondSim = getInitialConditions(data);
% plotFitResultSim(data, modelParamsSim, initCondSim, "cov", "reduced", errBins)
% 
% for i = 1:numel(modelParams)
%     fprintf("Data: %.4f, Sim: %4f \n", modelParams(i), modelParamsSim(i))
% end
% 
% % NLL for each trial
% nllsSimData = getNLLForEachTrial(data, initCondSim, "cov", "reduced", modelParams);
% 
% figure
% hold on
% BinEdges = 38:0.5:50;
% histogram(nllsData, BinEdges, DisplayName="Data")
% histogram(nllsSimData, BinEdges, DisplayName="Simulated Data")
% xlabel("NLL")
% ylabel("Count")
% title("NLL for each trial")
% legend

% %% Plot fit results (Full)
% errBins = -90:3:90;
% [~, idx] = min(dataToSave.resultFullCov.f);
% modelParams = dataToSave.resultFullCov.x(idx, :);
% plotFitResult(formattedData, modelParams, initCond, "cov", "full", errBins)
% 
% % Simulation - GR too high (red flag)
% % There is something about guess rate that impacts recoverability
% [~, idx] = min(dataToSave.resultFullCovSim.f);
% modelParams = dataToSave.resultFullCovSim.x(idx, :);
% data = GenerateCovData( modelParams, n_uncertainty_levels, 'full', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3) );
% initCondSim = getInitialConditions(data);
% plotFitResultSim(data, modelParams, initCondSim, "cov", "full", errBins)

% %% Plot fit results (Full)
% errBins = -90:3:90;
% [~, idx] = min(dataToSave.resultFullJumboCov.f);
% modelParams = dataToSave.resultFullJumboCov.x(idx, :);
% plotFitResult(formattedData, modelParams, initCond, "cov", "full_jumbo", errBins)
% 
% % Simulation - GR too high (red flag)
% % There is something about guess rate that impacts recoverability
% [~, idx] = min(dataToSave.resultFullJumboCovSim.f);
% modelParams = dataToSave.resultFullJumboCovSim.x(idx, :);
% data = GenerateCovData( modelParams, n_uncertainty_levels, 'full_jumbo', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3) );
% initCondSim = getInitialConditions(data);
% plotFitResultSim(data, modelParams, initCondSim, "cov", "full_jumbo", errBins)


%% Plot fit results (Reduced)
map.ind.reduced = "resultReducedInd";
map.ind.full    = "resultFullInd";
map.ind.jumbo   = "resultFullJumboInd";

map.cov.reduced = "resultReducedCov";
map.cov.full    = "resultFullCov";
map.cov.jumbo   = "resultFullJumboCov";

modelType = "ind";
fitType = "reduced";
key = map.(modelType).(fitType);
keySim = key + "Sim";

errBins = -90:3:90;
[~, idx] = min(dataToSave.(key).f);
modelParams = dataToSave.(key).x(idx, :);
% plotFitResult(formattedData, modelParams, initCond, modelType, fitType, errBins)
% NLL for each trial
nllsData = getNLLForEachTrial(formattedData, initCond, modelType, fitType, modelParams);

% Simulation - GR too high (red flag)
% There is something about guess rate that impacts recoverability
[~, idx] = min(dataToSave.(keySim).f);
modelParamsSim = dataToSave.(keySim).x(idx, :);
if modelType == "ind"
    data = GenerateIndData( modelParamsSim, n_uncertainty_levels, fitType, ...
        formattedData.orientations', size(formattedData.theta_resp_all, 3) );
elseif modelType == "cov"
    data = GenerateCovData( modelParamsSim, n_uncertainty_levels, fitType, ...
        formattedData.orientations', size(formattedData.theta_resp_all, 3) );
end

initCondSim = getInitialConditions(data);
% plotFitResultSim(data, modelParamsSim, initCondSim, modelType, fitType, errBins)

for i = 1:numel(modelParams)
    fprintf("Data: %.4f, Sim: %4f \n", modelParams(i), modelParamsSim(i))
end

% NLL for each trial
nllsSimData = getNLLForEachTrial(data, initCondSim, modelType, fitType, modelParams);

figure
hold on
BinEdges = 38:0.5:50;
histogram( nllsSimData - nllsData, DisplayName="del NLL(Data) - NLL(Sim)")
% histogram(nllsData, BinEdges, DisplayName="Data")
% histogram(nllsSimData, BinEdges, DisplayName="Simulated Data")
xline(0, LineStyle="--", LineWidth=1, HandleVisibility="off")
xlabel("del NLL(Data) - NLL(Sim)")
ylabel("Count")
title("del NLL(Data) - NLL(Sim) for each trial")
legend


%% Ind vs Cov
modelType = "ind";
fitType = "full";
key = map.(modelType).(fitType);

[~, idx] = min(dataToSave.(key).f);
modelParams = dataToSave.(key).x(idx, :);
% NLL for each trial
nllsDataInd = getNLLForEachTrial(formattedData, initCond, modelType, fitType, modelParams);

modelType = "cov";
fitType = "full";
key = map.(modelType).(fitType);

[~, idx] = min(dataToSave.(key).f);
modelParams = dataToSave.(key).x(idx, :);
% NLL for each trial
nllsDataCov = getNLLForEachTrial(formattedData, initCond, modelType, fitType, modelParams);

figure
hold on
BinEdges = 38:0.5:50;
histogram( nllsDataCov - nllsDataInd, DisplayName="del NLL(Cov) - NLL(Ind)")
% histogram(nllsData, BinEdges, DisplayName="Data")
% histogram(nllsSimData, BinEdges, DisplayName="Simulated Data")
xline(0, LineStyle="--", LineWidth=1, HandleVisibility="off")
xlabel("del NLL(Cov) - NLL(Ind)")
ylabel("Count")
title("del NLL(Cov) - NLL(Ind) for each trial")
legend