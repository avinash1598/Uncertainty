close all
clear all
restoredefaultpath

% ITR 191
warning("Make sure generative model is correct. While generating simulated " + ...
    "data guess rate should only be used for low confidence reports.")

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsNoBin/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsNoBin')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

dataFileName  = ["COR31.mat", "COR33.mat", "CORNFB01.mat", "CORNFB02.mat", "CORNFB03.mat"];
subjects      = ["Tien", "Akash", "Yichao", "Jonathan", "David"];
modelTypes    = ["cov", "ind", "singlyStochastic"];
fltSessionIdx = [0, 0, 0, 2, 1]; % Last one is for David - 1

%%
fitResults = {};

for i = 1:numel(subjects)
    fprintf("Finding parameters for subject: %s", subjects(i))

    % Load and format data
    clear expData
    clear formattedData
    clear initCond

    expData = load('../Data/'+dataFileName(i));

    fltData       = expData.dat( expData.dat.session > fltSessionIdx(i) , :);  % TODO: change session number
    f.dat         = fltData;
    formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

    initCond             = getInitialConditions(formattedData);
    n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now

    % Optimize
    for mIdx = 1:numel(modelTypes)
        fprintf("Optimizing Subject: %s, Model Type: %s", subjects(i), modelTypes(mIdx))

        modelType = modelTypes(mIdx);
        optParams.nStarts = 30;
        result = Optimize(formattedData, initCond, modelType, [], optParams, "full");

        [nll,idx] = min(result.f);
        fitParams = result.x(idx, :);
        
        fitResults.(subjects(i)).(modelType).nll    = nll;
        fitResults.(subjects(i)).(modelType).params = fitParams;

        fprintf("Subject: %s, Model Type: %s, NLL: %.2f", ...
            subjects(i), modelType, nll)
    end

    save("fitResultsAllSubjects.mat", "fitResults");
end

%% TODO: compute AIC/BIC after this and then save it to the same file.
load("fitResultsAllSubjects.mat");

covNLLs = zeros(1, numel(subjects));
indNLLs = zeros(1, numel(subjects));
singlyStochasticNLLs = zeros(1, numel(subjects));

for i = 1:numel(subjects)
    covNLLs(i) = fitResults.(subjects(i)).cov.nll;
    indNLLs(i) = fitResults.(subjects(i)).ind.nll;
    singlyStochasticNLLs(i) = fitResults.(subjects(i)).singlyStochastic.nll;
    
    % COV - SS
    k1 = 13; k2 = 11; n = 1728;
    [aic_cov_ss, bic_cov_ss]  = computeAIC_BIC(covNLLs(i), singlyStochasticNLLs(i), k1, k2, n); 
    fitResults.(subjects(i)).aic_cov_ss = aic_cov_ss;
    fitResults.(subjects(i)).bic_cov_ss = bic_cov_ss;

    % COV - IND
    k1 = 13; k2 = 13; n = 1728;
    [aic_cov_ind, bic_cov_ind]  = computeAIC_BIC(covNLLs(i), indNLLs(i), k1, k2, n); 
    fitResults.(subjects(i)).aic_cov_ind = aic_cov_ind;
    fitResults.(subjects(i)).bic_cov_ind = bic_cov_ind;

    % IND - SS
    k1 = 13; k2 = 11; n = 1728;
    [aic_ind_ss, bic_ind_ss]  = computeAIC_BIC(indNLLs(i), singlyStochasticNLLs(i), k1, k2, n); 
    fitResults.(subjects(i)).aic_ind_ss = aic_ind_ss;
    fitResults.(subjects(i)).bic_ind_ss = bic_ind_ss;

end

save("fitResultsAllSubjects.mat", "fitResults");

function [delAIC, delBIC] = computeAIC_BIC(nll1, nll2, k1, k2, n)

% AIC
AIC_1 = 2*k1 + 2*nll1;
AIC_2 = 2*k2 + 2*nll2;

% BIC
BIC_1 = k1*log(n) + 2*nll1;
BIC_2 = k2*log(n) + 2*nll2;

delAIC = AIC_1 - AIC_2;
delBIC = BIC_1 - BIC_2;

end
