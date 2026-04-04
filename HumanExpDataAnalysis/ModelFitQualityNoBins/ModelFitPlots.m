close all
clear all
restoredefaultpath

% GR: 0.05

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


bestCaseNLLs = [8126.6299, 7898.2796, 7883.5219]; %  8073.6588, 
worstCaseNLLs = [13727.4158, 14280.8876, 10980.2732]; %   12789.1775,  in actual this is a distribution
worstCaseNLLs2 = 10171.1878; % same for all
subjects = ["Akash", "Jonathan", "Yichao"]; % "Akash", "Tien", 

red   = [1 0 0];
blue  = [0 0 1];
green = [0 0.6 0];

% lighter shades (mix with white)
light = @(c,a) c + (1-c)*a; % a in [0,1]

red_light   = light(red, 0.5);
blue_light  = light(blue, 0.5);
green_light = light(green, 0.5);

for i = 1:numel(subjects)
    figure
    load(subjects(i) + "_cov_fit.mat");
    nllsCov = sort( resultFullCov.f, 1);
    
    load(subjects(i) + "_ind_fit.mat");
    nllsInd = sort( resultFullCov.f, 1);

    load(subjects(i) + "_singlyStochastic_fit.mat");
    nllsSS = sort( resultFullCov.f, 1);

    % COV - SS
    k1 = 13; k2 = 11; n = 1728;
    [aic_cov_ss, bic_cov_ss]  = computeAIC_BIC(nllsCov, nllsSS, k1, k2, n); 

    % COV - IND
    k1 = 13; k2 = 13; n = 1728;
    [aic_cov_ind, bic_cov_ind]  = computeAIC_BIC(nllsCov, nllsInd, k1, k2, n); 

    % IND - SS
    k1 = 13; k2 = 11; n = 1728;
    [aic_ind_ss, bic_ind_ss]  = computeAIC_BIC(nllsInd, nllsSS, k1, k2, n); 

    subplot(2, 1, 1)
    hold on
    histogram(nllsCov(1:15), BinWidth=2, DisplayName="Cov")
    histogram(nllsInd(1:15), BinWidth=2, DisplayName="Ind")
    histogram(nllsSS(1:15), BinWidth=2, DisplayName="Singly Stochastic")

    xline(bestCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Best case NLL", Color="green")
    xline(worstCaseNLLs(i), LineWidth=1.5, LineStyle="--", DisplayName="Worst case NLL", Color="red")

    xlabel("NLLs")
    ylabel("count")
    title("Guess rate bound: 1")
    legend
    ylim([0 15])
    xlim([bestCaseNLLs(i)-100, worstCaseNLLs(i) + 100])
    hold off
    
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

    subplot(2, 1, 2)
    hold on
    
    histogram(aic_cov_ss(1:15), 'BinWidth',2, FaceColor=red, DisplayName="AIC: COV- SS")
    histogram(aic_cov_ind(1:15), 'BinWidth',2, FaceColor=green, DisplayName="AIC: COV - IND")
    histogram(aic_ind_ss(1:15), 'BinWidth',2, FaceColor=blue, DisplayName="AIC: IND - SS")
    
    histogram(bic_cov_ss(1:15), 'BinWidth',2, 'FaceColor',red_light, DisplayName="BIC: COV - SS")
    histogram(bic_cov_ind(1:15), 'BinWidth',2, 'FaceColor',green_light, DisplayName="BIC: COV - IND")
    histogram(bic_ind_ss(1:15), 'BinWidth',2, 'FaceColor',blue_light, DisplayName="BIC: IND - SS")
    
    xlabel("AIC/BIC diff")
    ylabel("count")
    title("Diff (DS - SS)")
    legend
    ylim([0 15])
    xlim([-400 100])
    hold off

    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

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
