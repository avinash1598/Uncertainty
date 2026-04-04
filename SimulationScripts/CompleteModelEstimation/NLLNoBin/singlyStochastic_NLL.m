close all
clear all
restoredefaultpath

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
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsNoBin/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsNoBin')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% n_uncertainty_levels = 6;
% orientations     = 0:15:175; 
% ori_scale        = 0.67;
% b                = linspace(0.1, 1.5, 6); 
% a                = ori_scale.*b; 
% biasAmp          = 0.5; 
% shape            = 2;
% scale            = 0.5; %0.5;
% sigma_meta       = 0.2;
% Cc               = 0.5; 
% guessRate        = 0;
                             
orientations     = 0:15:175; 
n_uncertainty_levels = 6;
ori_scale        = 0.1626; %0.23;
b                = [6.6680 9.8866 13.4934 21.0658 33.0866 37.9527];
a                = ori_scale.*b;
biasAmp          = 1.5791; %10; %0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
sigma_meta       = 31.1465;
Cc               = 0.0503;   
guessRate        = 0.3; %0.0350;

modelParams = [b sigma_meta Cc guessRate ori_scale biasAmp];

%VV Imp: Remember to simulate guess rate only for low confidence trials
data  = GenerateData( ...
    modelParams, n_uncertainty_levels, 'full', orientations, 24);

initCond = getInitialConditions(data);

% errBins = -90:2:90;
% plotFitResultSim(data, modelParams, initCond, "singlyStochastic", "full", errBins)

%% Optimize
optParams.nStarts = 1;
result = Optimize(data, initCond, "singlyStochastic", [], optParams, "full");

%%
[~,idx] = min(result.f);
fitParams = result.x(idx, :);

for i=1:numel(result.x(idx, :))
    fprintf("GT: %.4f, Fit: %.4f \n", modelParams(i), fitParams(i))
end

% NLL computed using binned method
trlNLLs = getNLLForEachTrial(data, initCond, "singlyStochastic", "full", fitParams, []);
nll_binned_method = trlNLLs;

fprintf("\n\nNLL (fit) %.4f, NLL (binned method): %.4f) \n\n", result.f(idx), sum(nll_binned_method));

%% Plot (Simulated data and then optimization fit)
errBins = -90:2:90;
% plotFitResultSim(data, modelParams, initCond, "cov", "full", errBins)
plotFitResultSim(data, fitParams, initCond, "singlyStochastic", "full", errBins)


