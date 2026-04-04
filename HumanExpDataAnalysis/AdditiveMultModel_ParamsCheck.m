% TODO: do ti for all subjects

close all
clear all
restoredefaultpath

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
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsModel_v4/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScripts_v4/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsModel_v3/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/OptimizationScripts_v3/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/SimulationScripts/CompleteModelEstimation/GenerateSimulatedData')

% addpath('C:\Users\avinash1598\Desktop\Uncertainty\HumanExpDataAnalysis\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\LLScriptsUtils\LLScriptsTrialData\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\PlotUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\Utils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\OptimizationUtils\')
% addpath('C:\Users\avinash1598\Desktop\Uncertainty\SimulationScripts\CompleteModelEstimation\GenerateSimulatedData')

expData            = load('./Data/COR33.mat');     % Akash
% expData            = load('./Data/COR31.mat');   % Tien
% expData            = load('./Data/COR32.mat');   % Jiaming
% expData            = load('./Data/CORNFB02.mat');% Jonathan
% expData            = load('./Data/CORNFB01.mat');  % Yichao
% For subject 6 => David: make sure to change the orientation dependent
% error

fltSessionIdx = [0, 0, 0, 0, 2, 1];
fltData       = expData.dat( expData.dat.session > 0 , :);  % TODO: change session number
f.dat         = fltData;
formattedData = formatExpData(f, false, false); % no de-baising, work with raw errors

initCond      = getInitialConditions(formattedData);

n_uncertainty_levels = formattedData.n_uncertainty_levels; % Hard code for now

%%
% load('./Data/ModelFitQualityTest_Yichao_v4.mat');
load('./Data/ModelFitQualityTest_Akash_v4.mat');
% load('./Data/ModelFitQualityTest_Tien_v4.mat');
% load('./Data/ModelFitQualityTest_Jonathan_v4.mat');
errBins = -90:3:90;
[~, idx] = min(dataToSave.result.f);
modelParams = dataToSave.result.x(idx, :);
plotFitResult_v4(formattedData, modelParams, initCond, errBins)

opt_param_sigma_ext         = dataToSave.result.x(idx, 1:n_uncertainty_levels);
opt_param_sigma_b           = dataToSave.result.x(idx, n_uncertainty_levels+1:2*n_uncertainty_levels);
opt_param_scale             = dataToSave.result.x(idx, 2*n_uncertainty_levels + 1);
opt_param_sigma_meta        = dataToSave.result.x(idx, 2*n_uncertainty_levels + 2);
opt_param_Cc                = dataToSave.result.x(idx, 2*n_uncertainty_levels + 3);
opt_param_guessrate         = dataToSave.result.x(idx, 2*n_uncertainty_levels + 4);
opt_param_ori_scale         = dataToSave.result.x(idx, 2*n_uncertainty_levels + 5);
opt_param_bias              = dataToSave.result.x(idx, 2*n_uncertainty_levels + 6);

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("sigma_ext: Fit: %.4f C: %.3f, S: %.3f, D: %.3f\n", ...
        opt_param_sigma_ext(i), ...
        formattedData.uncertaintyVals(i, 1), ...
        formattedData.uncertaintyVals(i, 2), ...
        formattedData.uncertaintyVals(i, 3))
end

T1 = table( ...
    opt_param_sigma_ext(:), ...
    formattedData.uncertaintyVals(:,1), ...
    formattedData.uncertaintyVals(:,2), ...
    formattedData.uncertaintyVals(:,3), ...
    'VariableNames', {'sigma_ext_fit','contrast','spread','duration'} ...
);

% Mean internal noise
sigma_int_stim  = opt_param_sigma_b' + (opt_param_sigma_b'.*opt_param_ori_scale).*( ...
    abs(sind(initCond.sigma_m_shape1*formattedData.orientations' - initCond.sigma_m_shape2)));
sigma_int_stim_reduced = sqrt( mean( sigma_int_stim.^2, 2) );

% shape*scale = sigma_int^2
% *opt_param_scale /sqrt(formattedData.stimulusEnergy(i))
% opt_param_sigma_b(i)
for i =1:n_uncertainty_levels
    fprintf("b (sigma_int): Fit: %.4f C: %.3f, S: %.3f, D: %.3f\n", ...
        sigma_int_stim_reduced(i), ...
        formattedData.uncertaintyVals(i, 1), ...
        formattedData.uncertaintyVals(i, 2), ...
        formattedData.uncertaintyVals(i, 3))
end

T2 = table( ...
    sigma_int_stim_reduced(:), ...
    formattedData.uncertaintyVals(:,1), ...
    formattedData.uncertaintyVals(:,2), ...
    formattedData.uncertaintyVals(:,3), ...
    'VariableNames', {'sigma_int_fit','contrast','spread','duration'} ...
);

fprintf("Scale Fit: %.4f \n", opt_param_scale)
fprintf("Meta Fit: %.4f \n", opt_param_sigma_meta)
fprintf("Cc Fit: %.4f \n", opt_param_Cc)
fprintf("GR Fit: %.4f \n", opt_param_guessrate)
fprintf("Ori scale Fit: %.4f \n", opt_param_ori_scale)
fprintf("Bias Amp Fit: %.4f \n", opt_param_bias)


%% Plot params
uniqContrasts = unique( f.dat.stimContrast );
uniqSpreads   = unique( f.dat.stimSpread );
uniqDurations = unique( f.dat.stimDur );

pltIdx = 1;
figure

% Spread manipulation
for c=1:numel(uniqContrasts)
    for d=1:numel(uniqDurations)

        fltT1 = T1(T1.contrast == uniqContrasts(c) & T1.duration == uniqDurations(d), :); % external noise
        fltT2 = T2(T2.contrast == uniqContrasts(c) & T2.duration == uniqDurations(d), :); % internal noise

        if ~isempty(fltT1)
            subplot(3, 4, pltIdx)
            plot(fltT1.spread, fltT1.sigma_ext_fit)
            xlabel("Spread")
            ylabel("\sigma_ext")
            ylim([0 70])
            title(sprintf("C:%.3f, D:%.3f", uniqContrasts(c), uniqDurations(d)))
    
            subplot(3, 4, 4+pltIdx)
            plot(fltT2.spread, fltT2.sigma_int_fit)
            xlabel("Spread")
            ylabel("\sigma_int")
            ylim([0 70])
    
            pltIdx = pltIdx + 1;
        end

    end
end


pltIdx = 1;
figure

% Contrast manipulation
for s=1:numel(uniqSpreads)
    for d=1:numel(uniqDurations)

        fltT1 = T1(T1.spread == uniqSpreads(s) & T1.duration == uniqDurations(d), :); % external noise
        fltT2 = T2(T2.spread == uniqSpreads(s) & T2.duration == uniqDurations(d), :); % internal noise
        
        if size(fltT1, 1) > 1
            subplot(3, 4, pltIdx)
            plot(fltT1.contrast, fltT1.sigma_ext_fit)
            xlabel("Contrast")
            ylabel("\sigma_ext")
            ylim([0 70])
            title(sprintf("S:%.3f, D:%.3f", uniqSpreads(s), uniqDurations(d)))
            
            subplot(3, 4, 4+pltIdx)
            plot(fltT2.contrast, fltT2.sigma_int_fit)
            xlabel("Contrast")
            ylabel("\sigma_int")
            ylim([0 70])
    
            pltIdx = pltIdx + 1;
        end

    end
end


pltIdx = 1;
figure

% Duration manipulation
for s=1:numel(uniqSpreads)
    for c=1:numel(uniqContrasts)

        fltT1 = T1(T1.spread == uniqSpreads(s) & T1.contrast == uniqContrasts(c), :); % external noise
        fltT2 = T2(T2.spread == uniqSpreads(s) & T2.contrast == uniqContrasts(c), :); % internal noise
        
        if size(fltT1, 1) > 1
            subplot(3, 4, pltIdx)
            plot(fltT1.duration, fltT1.sigma_ext_fit)
            xlabel("Duration")
            ylabel("\sigma_ext")
            ylim([0 70])
            title(sprintf("S:%.3f, C:%.3f", uniqSpreads(s), uniqContrasts(d)))
            
            subplot(3, 4, 4+pltIdx)
            plot(fltT2.duration, fltT2.sigma_int_fit)
            xlabel("Duration")
            ylabel("\sigma_int")
            ylim([0 70])
    
            pltIdx = pltIdx + 1;
        end

    end
end
