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

load("fitResultsAllSubjects.mat");

subjects      = ["Tien", "Akash", "Yichao", "Jonathan"];
modelTypes    = ["cov", "ind", "singlyStochastic"];

covNLLs = zeros(1, numel(subjects));
indNLLs = zeros(1, numel(subjects));
singlyStochastic = zeros(1, numel(subjects));

aic_cov_ss_list = zeros(1, numel(subjects));
bic_cov_ss_list = zeros(1, numel(subjects));
aic_cov_ind = zeros(1, numel(subjects));
bic_cov_ind = zeros(1, numel(subjects));
aic_ind_ss  = zeros(1, numel(subjects));
bic_ind_ss  = zeros(1, numel(subjects));

for i = 1:numel(subjects)
    covNLLs(i) = fitResults.(subjects(i)).cov.nll;
    indNLLs(i) = fitResults.(subjects(i)).ind.nll;
    singlyStochastic(i) = fitResults.(subjects(i)).singlyStochastic.nll;

    aic_cov_ss_list(i) = fitResults.(subjects(i)).aic_cov_ss;
    bic_cov_ss_list(i) = fitResults.(subjects(i)).bic_cov_ss;
    aic_cov_ind(i)     = fitResults.(subjects(i)).aic_cov_ind;
    bic_cov_ind(i)     = fitResults.(subjects(i)).bic_cov_ind;
    aic_ind_ss(i)      = fitResults.(subjects(i)).aic_ind_ss;
    bic_ind_ss(i)      = fitResults.(subjects(i)).bic_ind_ss;

end

figure

subplot(2, 2, 1)
x = 1:numel(subjects); % or your actual x values
Y = [covNLLs(:), indNLLs(:), singlyStochastic(:)]; % combine columns
bar(x, Y)
legend({'cov','ind','SS'})
xlabel("Subjects")
ylabel("NLL")
ylim([8400 max(Y(:))+100])

% Is Ind and Cov better than SS?
subplot(2, 2, 2)
x = 1:numel(subjects); % or your actual x values
Y = [aic_cov_ss_list(:), bic_cov_ss_list(:)]; % combine columns
bar(x, Y)
legend({'AIC','BIC'})
xlabel("Subjects")
ylabel("diff (COV - SS)")
title("Is cov better than SS?")

% Is Ind better than SS?
subplot(2, 2, 3)
x = 1:numel(subjects); % or your actual x values
Y = [aic_ind_ss(:), bic_ind_ss(:)]; % combine columns
bar(x, Y)
legend({'AIC','BIC'})
xlabel("Subjects")
ylabel("diff (IND - SS)")
title("Is ind better than SS?")

% Is Cov better than Ind?
subplot(2, 2, 4)
x = 1:numel(subjects); % or your actual x values
Y = [aic_cov_ind(:), bic_cov_ind(:)]; % combine columns
bar(x, Y)
legend({'AIC','BIC'})
xlabel("Subjects")
ylabel("diff (COV - IND)")
title("Is cov better than ind?")


%% Plot params and fit

% visualize all parameters

% cov
paramsName = ["b L1",...
    "b L2", ...
    "b L3", ...
    "b L4", ...
    "b L5", ...
    "b L6", ...
    "scale", "sigma_meta", "Cc", "exp", "GR", "Ori scale", "biasAmp"];

% ind
paramsName = ["b L1",...
    "b L2", ...
    "b L3", ...
    "b L4", ...
    "b L5", ...
    "b L6", ...
    "mu_ln", "sigma_ln", "sigma_meta", "Cc", "GR", "Ori scale", "biasAmp"];

