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


bestCaseNLLs = [8126.6299, 8073.6588, 7898.2796, 7883.5219];
worstCaseNLLs = [13727.4158, 12789.1775, 14280.8876, 10980.2732]; %  in actual this is a distribution
worstCaseNLLs2 = 10171.1878; % same for all
subjects = ["Akash", "Tien", "Jonathan", "Yichao"];

for i = 1:numel(subjects)
    figure

    % No bound on guess rate
    load(subjects(i) + "_descriptive_model_fit.mat");
    doublyStochasticNLLList = dataToSave.covTotalNLLList;
    singlyStochasticNLLList = dataToSave.singlyStochasticTotalNLLList;
    
    itrNLL_DS = getNLLList(doublyStochasticNLLList);
    itrNLL_SS = getNLLList(singlyStochasticNLLList);

    diff1 = itrNLL_DS - itrNLL_SS;
    k1 = 42; k2 = 36; n = 1728;
    [aic1, bic1]  = computeAIC_BIC(itrNLL_DS, itrNLL_SS, k1, k2, n); 

    subplot(5, 1, 1)
    hold on
    histogram(itrNLL_DS(1:15), BinWidth=2, DisplayName="Doubly Stochastic")
    histogram(itrNLL_SS(1:15), BinWidth=2, DisplayName="Singly Stochastic")

    xline(bestCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Best case NLL", Color="green")
    xline(worstCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Worst case NLL", Color="red")

    xlabel("NLLs")
    ylabel("count")
    title("Guess rate bound: 1")
    legend
    ylim([0 15])
    hold off

    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

    load(subjects(i) + "_descriptive_model_fit_GR_0_15.mat");
    doublyStochasticNLLList = dataToSave.covTotalNLLList;
    singlyStochasticNLLList = dataToSave.singlyStochasticTotalNLLList;
    
    itrNLL_DS = getNLLList(doublyStochasticNLLList);
    itrNLL_SS = getNLLList(singlyStochasticNLLList);

    diff2 = itrNLL_DS - itrNLL_SS;
    k1 = 42; k2 = 36; n = 1728;
    [aic2, bic2]  = computeAIC_BIC(itrNLL_DS, itrNLL_SS, k1, k2, n); 

    subplot(5, 1, 2)
    hold on
    histogram(itrNLL_DS(1:15), BinWidth=2, DisplayName="Doubly Stochastic")
    histogram(itrNLL_SS(1:15), BinWidth=2, DisplayName="Singly Stochastic")

    xline(bestCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Best case NLL", Color="green")
    xline(worstCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Worst case NLL", Color="red")

    xlabel("NLLs")
    ylabel("count")
    title("Guess rate bound: 0.15")
    legend
    ylim([0 15])
    hold off

    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

    load(subjects(i) + "_descriptive_model_fit_GR_0_01.mat");
    doublyStochasticNLLList = dataToSave.covTotalNLLList;
    singlyStochasticNLLList = dataToSave.singlyStochasticTotalNLLList;
    
    itrNLL_DS = getNLLList(doublyStochasticNLLList);
    itrNLL_SS = getNLLList(singlyStochasticNLLList);

    diff3 = itrNLL_DS - itrNLL_SS;
    k1 = 42; k2 = 36; n = 1728;
    [aic3, bic3]  = computeAIC_BIC(itrNLL_DS, itrNLL_SS, k1, k2, n); 

    subplot(5, 1, 3)
    hold on
    histogram(itrNLL_DS(1:15), BinWidth=2, DisplayName="Doubly Stochastic")
    histogram(itrNLL_SS(1:15), BinWidth=2, DisplayName="Singly Stochastic")

    xline(bestCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Best case NLL", Color="green")
    xline(worstCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Worst case NLL", Color="red")

    xlabel("NLLs")
    ylabel("count")
    title("Guess rate bound: 0.01")
    legend
    ylim([0 15])
    hold off

    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';
    
    % NLL difference
    subplot(5, 1, 4)
    hold on
    histogram(diff1(1:15), BinWidth=2, DisplayName="GR bound: 1")
    histogram(diff2(1:15), BinWidth=2, DisplayName="GR bound: 0.15")
    histogram(diff3(1:15), BinWidth=2, DisplayName="GR bound: 0.01")

    xlabel("NLLs diff")
    ylabel("count")
    title("Diff (DS - SS)")
    legend
    ylim([0 15])
    hold off

    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

    % Param difference (AIC would probably be zero)
    subplot(5, 1, 5)
    hold on
    
    histogram(aic1(1:15), 'BinWidth',2, 'FaceColor',[1 0.4 0.4], DisplayName="AIC:  GR bound: 1")
    histogram(aic2(1:15), 'BinWidth',2, 'FaceColor',[1 0.2 0.2], DisplayName="AIC: GR bound: 0.15")
    histogram(aic3(1:15), 'BinWidth',2, 'FaceColor',[0.8 0 0], DisplayName="AIC: GR bound: 0.01")
    
    histogram(bic1(1:15), 'BinWidth',2, 'FaceColor',[0.4 0.4 1], DisplayName="BIC:  GR bound: 1")
    histogram(bic2(1:15), 'BinWidth',2, 'FaceColor',[0.2 0.2 1], DisplayName="BIC: GR bound: 0.15")
    histogram(bic3(1:15), 'BinWidth',2, 'FaceColor',[0 0 0.8], DisplayName="BIC: GR bound: 0.01")

    xlabel("AIC/BIC diff")
    ylabel("count")
    title("Diff (DS - SS)")
    legend
    ylim([0 15])
    hold off

    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

end


function nllList = getNLLList(nllList)
    sortedNLL = sort(nllList, 1); % sort each column ascending
    nllList = sum(sortedNLL, 2);
end

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

