close all
clear all

addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils/')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')

% modelData                            = load('../modelContOriData_cov.mat');
modelData                            = load('../modelContOriData.mat');
% modelData.data.resp_err_all          = modelData.data.resp_err_all;
% modelData.data.confidence_report_all = modelData.data.confidence_report_all;

%% Cross-validate
errBins   = -90:0.5:90;

optParams.nStarts = 30;
optParams.hyperParamC1 = 0;
optParams.hyperParamC2 = 0;
optParams.randomGuessModel = true;

cv_result = NLLCrossValidate(modelData.data, errBins, 5, 5, optParams, 'full');
save('./CV_Data/cross_validation_full_ind_model_fit_method_2.mat', 'cv_result'); % TODO: chnage it to ind model

%% NLL on test data
errBins   = -90:0.5:90;

cv_result = load('./CV_Data/cross_validation_full_ind_model_fit_method_2.mat');

optParams.nStarts = 20;
optParams.hyperParamC1 = 0;
optParams.hyperParamC2 = 0;
optParams.randomGuessModel = true;

nllData = computeNLL_CV(modelData.data, errBins, cv_result.cv_result, optParams, 'full'); 

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