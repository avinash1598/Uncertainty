close all
clear all

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils/')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')

% modelData                            = load('../modelContOriData_cov.mat');
modelData                            = load('../modelContOriData.mat');
modelData.data.orientations          = unique( modelData.data.stimOri(:) );
modelData.data.theta_true_all        = modelData.data.stimOri;
% modelData.data.resp_err_all          = modelData.data.err;
% modelData.data.confidence_report_all = modelData.data.confReport;
modelData.data.stdByOri              = squeeze( std(modelData.data.resp_err_all, 0, 3) );
modelData.data.madByOri              = squeeze( mad(modelData.data.resp_err_all, 1, 3) );

% modelData.data.resp_err_all          = modelData.data.err;
% modelData.data.confidence_report_all = modelData.data.confReport;

%% Cross-validate
errBins   = -90:0.5:90;

optParams.nStarts = 30;
optParams.hyperParamC1 = 0;
optParams.hyperParamC2 = 0;
optParams.randomGuessModel = true;

cv_result = NLLCrossValidate(modelData.data, errBins, 5, 5, optParams, 'full');
save('./CV_Data/cross_validation_full_ind_model_fit_method_2.mat', 'cv_result'); % TODO: chnage it to ind model

%% NLL on test data
errBins   = -90:0.1:90;

cv_result = load('./CV_Data/cross_validation_ind_model.mat');

fitType = "reduced";
optParams.nStarts = 30;
optParams.hyperParamC1 = 0;
optParams.hyperParamC2 = 0;
optParams.randomGuessModel = true;

nllData = computeNLL_CV(modelData.data, errBins, cv_result.cv_result, optParams, fitType); 

%
figure
subplot(2, 2, 1)
hold on
histogram(nllData.nllCovModel/100, DisplayName='cov model')
histogram(nllData.nllIndModel/100, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("Fit on test data")

subplot(2, 2, 2)
histogram(nllData.deltaNLL) 
ylabel("Count")
xlabel("delta NLL (cov - ind)")
legend
title("Fit on test data")

subplot(2, 2, 3)
hold on
bins = 0:30:300;
histogram(nllData.fvalsCov/100, DisplayName='cov model')
histogram(nllData.fvalsInd/100, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("fvals (train dataset)")

subplot(2, 2, 4)
hold on
bins = 0:30:300;
histogram(nllData.minfvalsCov/100, DisplayName='cov model')
histogram(nllData.minfvalsInd/100, DisplayName='ind model')
hold off
ylabel("Count")
xlabel("NLL")
legend
title("bestfit fvals (train dataset)")

%% Make plots
modelParams = nllData.paramsCovModel;
plotFitResult_guessrate(modelData.data, modelParams, "cov", errBins, false)

% Display parameters
for i =1:6
    fprintf("Fit: %.4f \n", modelParams(i))
end

fprintf("Scale Fit: %.4f \n", modelParams(i + 1))
fprintf("Meta Fit: %.4f \n", modelParams(i + 2))
fprintf("Cc Fit: %.4f \n", modelParams(i + 3))
fprintf("GR Fit: %.4f \n", modelParams(i + 4))
% fprintf("Ori scale Fit: %.4f \n", modelParams(i + 5))
% fprintf("Bias: %.4f \n", modelParams(i + 6))
disp('\n\n')

modelParams = nllData.paramsIndModel;
plotFitResult_guessrate(modelData.data, modelParams, "ind", errBins, false)

% Display parameters
for i =1:6
    fprintf("Fit: %.4f \n", modelParams(i))
end

fprintf("Shape Fit: %.4f \n", modelParams(i + 1))
fprintf("Scale Fit: %.4f \n", modelParams(i + 2))
fprintf("Meta Fit: %.4f \n", modelParams(i + 3))
fprintf("Cc Fit: %.4f \n", modelParams(i + 4))
fprintf("GR Fit: %.4f \n", modelParams(i + 5))
% fprintf("Ori scale Fit: %.4f \n", modelParams(i + 6))
% fprintf("Bias: %.4f \n", modelParams(i + 7))
