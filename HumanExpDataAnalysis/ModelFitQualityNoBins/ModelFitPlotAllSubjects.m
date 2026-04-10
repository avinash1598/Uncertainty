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
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsNoBin/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScriptsNoBin')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

load("fitResultsAllSubjects.mat");

dataFileName  = ["COR31.mat", "COR33.mat", "CORNFB01.mat", "CORNFB02.mat", "CORNFB03.mat"]; % 
subjects      = ["Tien", "Akash", "Yichao", "Jonathan", "David"]; %
modelTypes    = ["cov", "ind", "singlyStochastic"];
fltSessionIdx = [0, 0, 0, 2, 1]; % Last one is for David - 1

covNLLs = zeros(1, numel(subjects));
indNLLs = zeros(1, numel(subjects));
singlyStochastic = zeros(1, numel(subjects));

aic_cov_ss_list = zeros(1, numel(subjects));
bic_cov_ss_list = zeros(1, numel(subjects));
aic_cov_ind = zeros(1, numel(subjects));
bic_cov_ind = zeros(1, numel(subjects));
aic_ind_ss  = zeros(1, numel(subjects));
bic_ind_ss  = zeros(1, numel(subjects));

paramsCov   = zeros(numel(subjects), 13);
paramsInd   = zeros(numel(subjects), 13);

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

    paramsCov(i, :)    = fitResults.(subjects(i)).cov.params;
    paramsInd(i, :)    = fitResults.(subjects(i)).ind.params;
end

figure

subplot(2, 2, 1)
x = 1:numel(subjects); % or your actual x values
Y = [covNLLs(:), indNLLs(:), singlyStochastic(:)]; % combine columns
bar(x, Y)
legend({'cov','ind','SS'})
xlabel("Subjects")
ylabel("NLL")
ylim([8000 max(Y(:))+100])

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

% % cov
paramsCovName = ["b1",...
    "b2", ...
    "b3", ...
    "b4", ...
    "b5", ...
    "b6", ...
    "scale", "sigma meta", "Cc", "exponent", "Guess Rate", "Ori Scale", "biasAmp"];

% ind
paramsIndName = ["b1",...
    "b2", ...
    "b3", ...
    "b4", ...
    "b5", ...
    "b6", ...
    "mu ln", "var ln", "sigma meta", "Cc", "Guess Rate", "Ori scale", "biasAmp"];

figure

for i=7:13
    subplot(3, 7, i - 7 + 1)
    x = 1:numel(subjects); % or your actual x values
    Y = paramsCov(:,i); % combine columns
    bar(x, Y)
    xlabel("Subjects")
    title(paramsCovName(i))
end

for i=7:13
    subplot(3, 7, 7 + i - 7 + 1)
    x = 1:numel(subjects); % or your actual x values
    Y = paramsInd(:,i); % combine columns
    bar(x, Y)
    xlabel("Subjects")
    title(paramsIndName(i))
end

annotation('textbox',[0 0.94 1 0.05],'String','Cov (top), Ind (bottom)','EdgeColor','none','HorizontalAlignment','center')

%% Plot visulization of parameters
covScaleParamIdx = 7;
covSigmaMetaParamIdx = 8;
covCcParamIdx = 9;
covExponentParamIdx = 10;
covGRParamIdx = 11;
covOriScaleParamIdx = 12;
covBiasAmpParamIdx = 13;

figure
for i = 1:numel(subjects)

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

    params = paramsCov(i,:);

    b = params(1);
    scale = params(covScaleParamIdx);
    exponent = params(covExponentParamIdx);
    oriScale = params(covOriScaleParamIdx);
    a = oriScale*b;
    stimOris = formattedData.orientations;

    sigma_m_shape1 = initCond.sigma_m_shape1;
    sigma_m_shape2 = initCond.sigma_m_shape2;
    sigma_s_stim = b + a.*(abs(sind(sigma_m_shape1*stimOris - sigma_m_shape2)));
    sigma_s_stim = max(sigma_s_stim); % Plot for the maximum sensory noise

    m = sigma_s_stim.^2; v = scale*(m^exponent);
    sigma2 = log(1 + v/m^2);  % variance
    mu = log(m) - 0.5*sigma2; %mean
    
    % transform to X (since Y = X^2)
    mu = mu/2;
    sigma2 = sigma2/4;
    
    x = linspace(0, sqrt(m)+3*sqrt( sqrt(v) ), 500);
    y = lognpdf(x, mu, sqrt(sigma2));
    
    subplot(3, 4, i)
    % plot + shaded area
    plot(x,y,'LineWidth',2, HandleVisibility='off'); hold on
    area(x,y,'FaceAlpha',0.3,'EdgeColor','none', HandleVisibility='off'); 
    xline(sigma_s_stim, LineStyle="--", LineWidth=1.5, DisplayName="\sigma_s"); hold off
    xlabel("\sigma_m")
    ylabel("P(\sigma_m)")
    legend
    title(sprintf("Doubly stochastic \n measurement noise"))

    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

    % Meta uncertainty
    m = sigma_s_stim;
    s = params(covSigmaMetaParamIdx); % sigma_meta
    sigma2 = log(1 + s^2 / m^2);
    mu = log(m) - 0.5*sigma2;

    x = linspace(0, m+4*s, 500);
    y = lognpdf(x, mu, sqrt(sigma2));
    
    subplot(3, 4, 4+i)
    % plot + shaded area
    plot(x,y,'LineWidth',2, HandleVisibility='off'); hold on
    area(x,y,'FaceAlpha',0.3,'EdgeColor','none', HandleVisibility='off'); 
    xline(m, LineStyle="--", LineWidth=1.5, DisplayName='\sigma_m'); hold off
    xlabel('\sigma_m est')
    ylabel("P(\sigma_m est)")
    legend
    title("Meta uncertainty")

    % Confidence variable
    Cc = params(covCcParamIdx);
    m = sigma_s_stim;
    s = params(covSigmaMetaParamIdx); % sigma_meta
    mu    = - log((m.^2) ./ sqrt(s.^2 + m.^2));
    sigma2 = (log((s.^2)./(m.^2) + 1));
    
    m = exp(mu + sigma2/2);
    s = sqrt( (exp(sigma2) - 1)*(exp(2*mu + sigma2)) );
    
    x = linspace(0, m+4*s, 500);
    y = lognpdf(x, mu, sqrt(sigma2));

    subplot(4, 4, 8+i)
    % plot + shaded area
    plot(x,y,'LineWidth',2, HandleVisibility='off'); hold on
    area(x,y,'FaceAlpha',0.3,'EdgeColor','none', HandleVisibility='off'); 
    xline(Cc, LineStyle="--", LineWidth=1.5, DisplayName='Cc'); hold off
    xlabel('V_c')
    ylabel("P(V_c)")
    legend
    title("Confidence variable")
end

%% Plot visulization of parameters
indMuLnParamIdx = 7;
indVarLnParamIdx = 8;
indSigmaMetaParamIdx = 9;
indCcParamIdx = 10;
indGRParamIdx = 11;
indOriScaleParamIdx = 12;
indBiasAmpParamIdx = 13;

figure
for i = 1:numel(subjects)

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

    params = paramsInd(i,:);
    
    b = params(1);
    muln = params(indMuLnParamIdx);
    sigmaln2 = params(indVarLnParamIdx);
    oriScale = params(indOriScaleParamIdx);
    a = oriScale*b;
    stimOris = formattedData.orientations;

    sigma_m_shape1 = initCond.sigma_m_shape1;
    sigma_m_shape2 = initCond.sigma_m_shape2;
    sigma_s_stim = b + a.*(abs(sind(sigma_m_shape1*stimOris - sigma_m_shape2)));
    sigma_s_stim = max(sigma_s_stim); % Plot for the maximum sensory noise
    
    m = muln + sigma_s_stim^2; v = sigmaln2; % desired mean & variance of lognormal distribution
    sigma2 = log(1 + v/m^2);  % variance
    mu = log(m) - 0.5*sigma2; %mean
    
    % transform to X (since Y = X^2)
    mu = mu/2;
    sigma2 = sigma2/4;
    
    x = linspace(0, sqrt(m)+2*sqrt( sqrt(v) ), 500);
    y = lognpdf(x, mu, sqrt(sigma2));
    
    subplot(4, 4, i)
    % plot + shaded area
    plot(x,y,'LineWidth',2, HandleVisibility='off'); hold on
    area(x,y,'FaceAlpha',0.3,'EdgeColor','none', HandleVisibility='off'); 
    xline(sigma_s_stim, LineStyle="--", LineWidth=1.5, DisplayName="\sigma_s"); hold off
    xlabel("\sigma_m")
    ylabel("P(\sigma_m)")
    legend
    title(sprintf("Doubly stochastic \n measurement noise"))

end

% %% Plot curve fits
% 
% figure
% 
% for i = 1:numel(subjects)
% 
%     fprintf("Finding parameters for subject: %s \n", subjects(i))
% 
%     % Load and format data
%     clear expData
%     clear formattedData
%     clear initCond
% 
%     expData = load('../Data/'+dataFileName(i));
% 
%     fltData       = expData.dat( expData.dat.session > fltSessionIdx(i) , :);  % TODO: change session number
%     f.dat         = fltData;
%     formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors
% 
%     initCond             = getInitialConditions(formattedData);
%     n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now
% 
%     % Cov fits only
%     modelParams = paramsCov(i, :);
%     
%     % Std Dev
%     pltIdx = [4, 4, 4*(i - 1) + 1];
%     plotMADs = false;
%     plotFitsMinimal(formattedData, modelParams, initCond, ...
%         "cov", "full", pltIdx, plotMADs)
% 
%     % MAD
%     pltIdx = [4, 4, 4*(i - 1) + 3];
%     plotMADs = true;
%     plotFitsMinimal(formattedData, modelParams, initCond, ...
%         "cov", "full", pltIdx, plotMADs)
% 
% end
