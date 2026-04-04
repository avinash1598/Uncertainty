close all
clear all
restoredefaultpath

% There is no ind or cov in discriptive model - only singlystochastic or
% doublystochastic

% ITR 191
warning("Make sure generative model is correct. While generating simulated " + ...
    "data guess rate should only be used for low confidence reports.")

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsNoBin\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\OptimizationScriptsNoBin\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsNoBin_dscrptv/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsNoBin_dscrptv')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% Bound GR: 1, 0.1, 0.15 0.01
% expData            = load('../Data/COR33.mat');    % Akash, -66.37, -311.35, -88, -153.6379
% expData            = load('../Data/COR31.mat');     % Tien -17.744, -57.416, -37,  -128.2702
% expData            = load('./Data/COR32.mat');     % Jiaming 
expData            = load('../Data/CORNFB02.mat');  % Jonathan -6.4672, -36.852, -20, -92.4352
% expData            = load('../Data/CORNFB01.mat');  % Yichao -33.263, -38, -34.981, -73.2772
% For subject 6 => David: make sure to change the orientation dependent
% error

fltSessionIdx = [0, 0, 0, 0, 2, 1]; % Last one is for David
sessionFltIdx = 2;

fltData       = expData.dat( expData.dat.session > sessionFltIdx, : );
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

% contrast, spread, duration
stimVals = formattedData.uncertaintyVals;

covParams              = zeros(6, 7);
singlyStochasticParams = zeros(6, 6);
covTotalNLL = 0;
singlyStochasticTotalNLL = 0;
covTotalNLLList = [];
singlyStochasticTotalNLLList = [];

for stimIdx = 1:size(stimVals, 1)

    fltData       = expData.dat( expData.dat.session > sessionFltIdx & ...
        (expData.dat.stimContrast == stimVals(stimIdx, 1)) & ...
        (expData.dat.stimSpread == stimVals(stimIdx, 2)) & ...
        (expData.dat.stimDur == stimVals(stimIdx, 3)), :);  % TODO: change session number
    
    f.dat         = fltData;
    formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors
    initCond      = getInitialConditions(formattedData);
    
    % ----------------------------------------
    optParams.nStarts = 30; %300;
    resultFullCov = Optimize(formattedData, initCond, "singlyStochastic", [], optParams, "full");

    [val,idx] = min(resultFullCov.f);
    fitParams = resultFullCov.x(idx, :);
    singlyStochasticTotalNLL = singlyStochasticTotalNLL + val;
    singlyStochasticParams(stimIdx, :) = fitParams;
    singlyStochasticTotalNLLList = [singlyStochasticTotalNLLList resultFullCov.f];

    optParams.nStarts = 30; %300;
    resultFullCov = Optimize(formattedData, initCond, "doublyStochastic", [], optParams, "full");
    
    [val,idx] = min(resultFullCov.f);
    fitParams = resultFullCov.x(idx, :);
    covTotalNLL = covTotalNLL + val;
    covParams(stimIdx, :) = fitParams;
    covTotalNLLList = [covTotalNLLList resultFullCov.f];

    % No need to test independent model here - Cov and Ind are literallly
    % same in discriptive sense
end

dataToSave.singlyStochasticTotalNLL = singlyStochasticTotalNLL;
dataToSave.singlyStochasticParams   = singlyStochasticParams;
dataToSave.singlyStochasticTotalNLLList = singlyStochasticTotalNLLList;
dataToSave.covTotalNLL              = covTotalNLL;
dataToSave.covParams                = covParams;
dataToSave.covTotalNLLList          = covTotalNLLList;

% save("Yichao_descriptive_model_fit_GR_0_01.mat", "dataToSave");

dataToSave.covTotalNLL - dataToSave.singlyStochasticTotalNLL

%% Plot parameters
% close all
% clear all

% load("Yichao_descriptive_model_fit_GR_0_15.mat", "dataToSave");
% covParams = dataToSave.covParams;
% singlyStochasticParams = dataToSave.singlyStochasticParams;

paramsNameCov = ["sigma_s", "scale", "sigma_meta", "Cc", "GR", "Ori scale", "biasAmp"];
paramsNameSinglyStochastic = ["sigma_s", "sigma_meta", "Cc", "GR", "Ori scale", "biasAmp"];

figure
for i=1:7
    subplot(2, 7, i)
    plot(1:6, covParams(:, i))
    title(paramsNameCov(i))

    if i <= 6
        subplot(2, 7, i+7)
        plot(1:6, singlyStochasticParams(:, i))
        title(paramsNameSinglyStochastic(i))
    end

end

z = covParams(:, 2)./covParams(:, 1);
x = covParams(:, 1);
[~, idx] = sort(x);
x = x(idx);
z = z(idx);

model = @(p,x) p(1)*x.^p(2); % + p(3);
p0 = [1 1];
p = lsqcurvefit(model, p0, x, z);

xfit = linspace(min(x), max(x), 100);
zfit = model(p, xfit);

figure; 
hold on
plot(x, z, 'o'); 
plot(xfit, zfit, '-'); hold off
% set(gca, 'YScale', 'log')


% %%
% [~,idx] = min(resultFullCov.f);
% fitParams = resultFullCov.x(idx, :);
% for i=1:numel(fitParams)
%     fprintf("Fit: %.4f \n", fitParams(i))
% end
% 
% 
% fprintf("\n\nNLL (fit) %.4f \n\n", resultFullCov.f(idx));

% Tien:

% NLL (fit) 1238.5245, 1400, 1440.7116, 1557, 1506.7884, 1555.1333
% COV: 8696 - absolute best with discriptive model
% SS: 1239.2272, 1400.7579, 1445.8390, 1557.9237, 1516.6439 , 1556.7034 = 8713

% Akash
% SS: 8909.6818
% IND: 8861.5184
% COV: 8851.8812 FULL: 8830.8654
% Best case: 8126.6299 

% Tien
% COV: 8987.6236 - FULL: 8767.9175, 8771.1925 
% IND: 8778.6861
% SS: 8786.5863
% Best: 8073

% Jonathan
% COV: 8767.3606 - 8773.2721 - 8766.8572 - 8763.0920 
% SS: 8772.8988 - 8786.0196
% IND: 8772.7061

% Yichao
% SS: 8601.538
% COV: 8571.7586
% COV FULL: 8565.6829

% GAMA to lognormal - try this
%%
% errBins = -90:2:90;
% plotFitResultSim(formattedData, fitParams, initCond, "cov", "full", errBins)

% % Simulated data
% [~, idx] = min(resultFullCov.f);
% params   = resultFullCov.x(idx, :);
% % params(end) = 0; % setting GR zero for now
% data     = GenerateCovData( params, n_uncertainty_levels, 'full', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3) );
% initCond = getInitialConditions(data);
% 
% optParams.nStarts = 300; %300
% resultFullCovSim = Optimize(data, initCond, "cov", [], optParams, 'full');

%% Ind model
% ----------------------------------------
% optParams.nStarts = 100; %300;
% resultFullInd = Optimize(formattedData, initCond, "ind", [], optParams, "full");

% % Simulated data
% [~, idx] = min(resultFullInd.f);
% params   = resultFullInd.x(idx, :);
% % params(end - 2) = 0;
% data     = GenerateIndData(params, n_uncertainty_levels, 'full', ...
%     formattedData.orientations', size(formattedData.theta_resp_all, 3));
% initCond = getInitialConditions(data);
% 
% optParams.nStarts = 1; %300
% resultFullIndSim = Optimize(data, initCond, "ind", [], optParams, 'full');


%% Singly Stochastic model
% optParams.nStarts = 1; %300;
% resultFullSinglyStochastic = Optimize(formattedData, initCond, "singlyStochastic", [], optParams, "full");

% SS: 8886.5879

% modelFitResults.resultFullCov              = resultFullCov;
% modelFitResults.resultFullInd              = resultFullInd;
% modelFitResults.resultFullSinglyStochastic = resultFullSinglyStochastic;
% 
% save('ModelFitQualityTest_Akash.mat', "modelFitResults");
