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
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsNoBin/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsNoBin')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')


addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

n_uncertainty_levels = 6;
orientations     = 0:15:175; 
ori_scale        = 0.67;
b                = linspace(1, 2.2, 6); 
a                = ori_scale.*b; 
biasAmp          = 0.5; 
scale            = 0.5; %0.5;
sigma_meta       = 0.6;
Cc               = 0.5; 
guessRate        = 0;
                             

% orientations     = 0:15:175; 
% n_uncertainty_levels = 6;
% ori_scale        = 0.21927; %0.23;
% b                = [10.793, 16.762, 17.282, 22.6, 29.52, 32.616 ];
% % b                = linspace(10, 35, n_uncertainty_levels); % 1.2 % Choose b such that average noise level ranges from low to high (relative to internal noise level)
% a                = ori_scale.*b;
% biasAmp          = 1.3903; %10; %0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% scale            = 196.05; %163.3809; %0.5;
% sigma_meta       = 11.559; %20; %9.7355; %0.6;
% Cc               = 0.046103; %0.0261; %0.5; 
% guessRate        = 0.3306; %0.33065; %0.3808;

% sigma_s_stim          = b' + (ori_scale.*b)'*(abs(sind(orientations - 90)));
% bias                  = biasAmp*sind(2*orientations); 
% sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
% sigma_s_reduced_model = [12.0081 18.7643 19.3554 25.0146 32.1753 35.2552];

% modelParams = [sigma_s_reduced_model scale sigma_meta Cc guessRate];
modelParams = [b scale sigma_meta Cc guessRate ori_scale biasAmp];
% modelParams = [b scale sigma_meta Cc guessRate a biasAmp];

%VV Imp: Remember to simulate guess rate only for low confidence trials
data  = GenerateCovData( ...
    modelParams, n_uncertainty_levels, 'full', orientations, 24);

initCond = getInitialConditions(data);

%% Optimize
optParams.nStarts = 1;
result = Optimize(data, initCond, "cov", [], optParams, "full");

%%
[~,idx] = min(result.f);
fitParams = result.x(idx, :);

for i=1:numel(result.x(idx, :))
    fprintf("GT: %.4f, Fit: %.4f \n", modelParams(i), fitParams(i))
end

% NLL computed using binned method
trlNLLs = getNLLForEachTrial(data, initCond, "cov", "full", fitParams, []);
nll_binned_method = trlNLLs;

fprintf("\n\nNLL (fit) %.4f, NLL (binned method): %.4f) \n\n", result.f(idx), sum(nll_binned_method));

%% Plot (Simulated data and then optimization fit)
errBins = -90:2:90;
plotFitResultSim(data, modelParams, initCond, "cov", "full", errBins)
plotFitResultSim(data, fitParams, initCond, "cov", "full", errBins)

