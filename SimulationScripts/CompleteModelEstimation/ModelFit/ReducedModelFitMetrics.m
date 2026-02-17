clear all
close all

addpath('data/')

% data = load('ModelFitMetrics.mat');
data = load('ModelFitMetrics_common_objective_t25.mat');
fitMetrics = data.dataToSave.fitInfo;
groundTruth = data.dataToSave.groundTruth;

tableSummary = groupsummary(fitMetrics, {'NumTrials'}, {'mean', 'std'}, {'paramsFitLoss', 'curveFitLoss'});

% Plot model fit losses

trials = tableSummary.NumTrials;
m_param_fit_loss = tableSummary.mean_paramsFitLoss;
std_param_fit_loss = tableSummary.std_paramsFitLoss./sqrt(tableSummary.GroupCount);
m_curve_fit_loss = tableSummary.mean_curveFitLoss;
std_curve_fit_loss = tableSummary.std_curveFitLoss;

figure
hold on
errorbar(1:numel(trials), m_param_fit_loss, std_param_fit_loss, DisplayName="Param fit", LineWidth=2)
errorbar(1:numel(trials), m_curve_fit_loss, std_curve_fit_loss, DisplayName="Curve fit", LineWidth=2)
xticks(1:numel(trials))
xticklabels(trials)
xlabel("Trials count")
ylabel("loss (mean)")
% ylim([0, 100])
legend
hold off

%% Plot parameters 
trlCnt = 250;
nParams = 12;
columnNames = fitMetrics.Properties.VariableNames;

figure

for i = 1:nParams
    x_ = fitMetrics{:, 1};
    fltIdx = (x_ == trlCnt);

    p_ = fitMetrics{:, i+2};
    p_ = p_(fltIdx);

    actual_param_val = groundTruth(i);
    subplot(4, 3, i)
    
    hold on
    histogram(p_, BinEdges=0:0.1:4)
    xline(actual_param_val, LineStyle="--", LineWidth=1.5)
    xlabel("Fit values")
    ylabel("count")
    title(columnNames{i+2})
    hold off
end

% %% Plot parameters 
% close all
% clear all
% 
% data = load('ModelFitMetrics_ntrials_2000.mat');
% fitMetrics = data.dataToSave.fitInfo;
% groundTruth = data.dataToSave.groundTruth;
% nParams = 12;
% 
% columnNames = fitMetrics.Properties.VariableNames;
% 
% figure
% 
% for i = 1:nParams
%     p_ = fitMetrics{:, i+2};
%     
%     actual_param_val = groundTruth(i);
%     subplot(4, 3, i)
%     
%     hold on
%     histogram(p_, BinEdges=0:0.1:5, DisplayName="Fit")
%     xline(actual_param_val, LineStyle="--", LineWidth=1.5, DisplayName="Actual val")
%     xlabel("Fit values")
%     ylabel("count")
%     title(columnNames{i+2})
%     legend
%     hold off
% end