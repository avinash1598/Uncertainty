restoredefaultpath
close all
clear all

% get(groot,"default")
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

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR31.mat'); % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR33.mat'); % Akash
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB01.mat'); % Yichao
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB02.mat');   % Jonathan
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB03.mat'); % David

subjectName = ["Tien", ...
    "Akash", ...
    "Jonathan", ...
     "Yichao"]; %" "Yichao"

dataFileNames = ["COR31.mat" "COR33.mat" "CORNFB01.mat" "CORNFB02.mat" ]; % "CORNFB03.mat"
fltSessionIdx = [0, 0, 0, 0, 2, 1];


figure

% Simulation vs data NLL comparison (Full Model)
for fIdx=1:numel(subjectName)

    % Add+SIM, Cov, Ind
    d1 = load("../Data/ModelFitQualityTest_" +  subjectName(fIdx) + "_v4.mat");
    
    % Cov
    d2 = load("../Data/ModelFitQualityTest_" +  subjectName(fIdx) + ".mat");
    
    fvalsAM  = d1.dataToSave.result.f;
    fvalsCov = d2.dataToSave.resultFullCov.f;
    fvalsInd = d2.dataToSave.resultFullInd.f;
    
    top100 = sort(fvalsAM, 'ascend'); top100 = top100(1:50); fvalsAM = top100;
    top100 = sort(fvalsCov, 'ascend'); top100 = top100(1:50); fvalsCov = top100;
    top100 = sort(fvalsInd, 'ascend'); top100 = top100(1:50); fvalsInd = top100;
    
    subplot(3, 4, fIdx)
    hold on
    histogram(fvalsAM, 'BinWidth', 2, DisplayName="ADD+MULT (full)"); %BinEdges,
    histogram(fvalsCov, 'BinWidth', 2, DisplayName="Cov");
    histogram(fvalsInd, 'BinWidth', 2, DisplayName="Ind");
    %xlabel("fvals (log)")
    xlabel("NLL")
    ylabel("count")
    title(sprintf("Model comparison (Subject: %d)", fIdx))
    legend
    hold off
    
    % p = ranksum(fvals, fvalsSim)
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';
    
    % AIC/BIC
    n = 1728;        % number of trials
    k_cov = size( d2.dataToSave.resultFullCov.x , 2); %12; % parameters (Cov model)
    k_ind = size( d2.dataToSave.resultFullInd.x , 2);      % parameters (Ind model)
    nll_cov = fvalsCov;
    nll_ind = fvalsInd;
    [delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);
    
    subplot(3,4,fIdx+4); hold on;
    histogram(delAIC, 'BinWidth', 2, 'DisplayName',' delta AIC');
    histogram(delBIC, 'BinWidth', 2, 'DisplayName','delta BIC');
    xlabel("delta (Cov - Ind)")
    title(sprintf("Model comparison (Subject: %d)", fIdx))
    legend;
    xline(0, LineStyle="--")
    xlim([-50 50])
    

    % AIC/BIC
    n = 1728;        % number of trials
    k_AM = size( d1.dataToSave.result.x , 2); %12; % parameters (Cov model)
    k_cov = size( d2.dataToSave.resultFullCov.x , 2);      % parameters (Ind model)
    nll_AM = fvalsAM;
    nll_cov = fvalsCov;
    [delAIC, delBIC] = computeAIC_BIC(nll_AM, nll_cov, k_AM, k_cov, n);
    
    subplot(3,4,fIdx+4+4); hold on;
    histogram(delAIC, 'BinWidth', 2, 'DisplayName',' delta AIC');
    histogram(delBIC, 'BinWidth', 2, 'DisplayName','delta BIC');
    xlabel("delta (Cov - Ind)")
    title(sprintf("Model comparison (Subject: %d)", fIdx))
    legend;
    xline(0, LineStyle="--")
    xlim([-150 150])
    
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
