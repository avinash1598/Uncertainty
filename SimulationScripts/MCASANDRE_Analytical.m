clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

% CASE1: sigma_e >> sigma_i : Low contrast  (7, 1)
% CASE2: sigma_e = sigma_i :                (1, 1)
% CASE3: sigma_e << sigma_i : High contrast (0.1, 1)

sigma_e_values = [0.1, 1, 5];
strArr = ["sigma_e << sigma_i", "sigma_e = sigma_i", "sigma_e >> sigma_i"];

for seIdx=1:numel(sigma_e_values)
    
    % Change in sigma_m
    stimVals = linspace(-20, 20, 100);
    sigma_m = [0.2, 7];
    
    modelParams.sigma_e                   = sigma_e_values(seIdx);
    modelParams.sigma_i                   = 1;
    modelParams.varGain                   = 0.5;
    modelParams.internalNoiseSamplesCnt   = 200;
    modelParams.Cd                        = 0;
    modelParams.Cc                        = 1;
    modelParams.sigma_m                   = 0.2;
    modelParams.sampleCnt                 = 100;
    
    figure
    sgtitle(strArr(seIdx))
    
    for idx=1:numel(sigma_m)
        modelParams.sigma_m = sigma_m(idx);
        choicePDFs = getLLhChoice_MCASANDRE(stimVals, modelParams);
        
        subplot(4, 2, 1)
        hold on
        psychFn = choicePDFs(3, :) + choicePDFs(4, :);
        plot(stimVals, psychFn, DisplayName="\sigma_m = " + sigma_m(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("% Choice 1")
        ylim([0, 1])
    
        hold off
        
        subplot(4, 2, 2)
        hold on
        
        prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
        plot(stimVals, prop_HC, DisplayName="\sigma_m = " + sigma_m(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("%HC")
        ylim([0, 1])
    
        hold off
    
    end
    
    
    % Change in Cc
    stimVals = linspace(-20, 20, 100);
    Cc = [1, 2];
    
    modelParams.sigma_e                   = sigma_e_values(seIdx);
    modelParams.sigma_i                   = 1;
    modelParams.varGain                   = 0.5;
    modelParams.internalNoiseSamplesCnt   = 200;
    modelParams.Cd                        = 0;
    modelParams.Cc                        = 1;
    modelParams.sigma_m                   = 0.2;
    modelParams.sampleCnt                 = 100;
    
    
    for idx=1:numel(Cc)
        modelParams.Cc = Cc(idx);
        choicePDFs = getLLhChoice_MCASANDRE(stimVals, modelParams);
        
        subplot(4, 2, 3)
        hold on
    
        psychFn = choicePDFs(3, :) + choicePDFs(4, :);
        plot(stimVals, psychFn, DisplayName="Cc = " + Cc(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("% Choice 1")
        ylim([0, 1])
    
        hold off
        
        subplot(4, 2, 4)
        hold on
        
        prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
        plot(stimVals, prop_HC, DisplayName="Cc = " + Cc(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("%HC")
        ylim([0, 1])
    
        hold off
    
    end
    
    
    % Change in sigma_d
    stimVals = linspace(-20, 20, 100);
    sigma_e = [0.1, 1, 3];
    
    modelParams.sigma_e                   = sigma_e_values(seIdx);
    modelParams.sigma_i                   = 1;
    modelParams.varGain                   = 0.5;
    modelParams.internalNoiseSamplesCnt   = 200;
    modelParams.Cd                        = 0;
    modelParams.Cc                        = 1;
    modelParams.sigma_m                   = 0.2;
    modelParams.sampleCnt                 = 100;
    
    
    for idx=1:numel(sigma_e)
        modelParams.sigma_e = sigma_e(idx);
        choicePDFs = getLLhChoice_MCASANDRE(stimVals, modelParams);
        
        subplot(4, 2, 5)
        hold on
        
        psychFn = choicePDFs(3, :) + choicePDFs(4, :);
        plot(stimVals, psychFn, DisplayName="\sigma_e = " + sigma_e(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("% Choice 1")
        ylim([0, 1])
    
        hold off
        
        subplot(4, 2, 6)
        hold on
        
        prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
        plot(stimVals, prop_HC, DisplayName="\sigma_e = " + sigma_e(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("%HC")
        ylim([0, 1])
    
        hold off
    
    end
    
    
    
    % Change in Cd
    stimVals = linspace(-20, 20, 100);
    Cd = [0, 1];
    
    modelParams.sigma_e                   = sigma_e_values(seIdx);
    modelParams.sigma_i                   = 1;
    modelParams.varGain                   = 0.5;
    modelParams.internalNoiseSamplesCnt   = 200;
    modelParams.Cd                        = 0;
    modelParams.Cc                        = 1;
    modelParams.sigma_m                   = 0.2;
    modelParams.sampleCnt                 = 100;
    
    
    for idx=1:numel(Cd)
        modelParams.Cd = Cd(idx);
        choicePDFs = getLLhChoice_MCASANDRE(stimVals, modelParams);
        
        subplot(4, 2, 7)
        hold on
        
        psychFn = choicePDFs(3, :) + choicePDFs(4, :);
        plot(stimVals, psychFn, DisplayName="Cd = " + Cd(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("% Choice 1")
        ylim([0, 1])
    
        hold off
        
        subplot(4, 2, 8)
        hold on
        
        prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
        plot(stimVals, prop_HC, DisplayName="Cd = " + Cd(idx), LineWidth=1.5, LineStyle="-")
        legend
        xlabel("Stimulus")
        ylabel("%HC")
        ylim([0, 1])
    
        hold off
    
    end

end


% % figure
% % subplot(2, 2, 1)
% 
% hold on
% 
% plot(stimVals, choicePDFs(1, :), DisplayName="D0, HC", LineWidth=1.5, LineStyle="-")
% plot(stimVals, choicePDFs(2, :), DisplayName="D0, LC", LineWidth=1.5, LineStyle="-")
% plot(stimVals, choicePDFs(3, :), DisplayName="D1, HC", LineWidth=1.5, LineStyle="-")
% plot(stimVals, choicePDFs(4, :), DisplayName="D1, LC", LineWidth=1.5, LineStyle="-")
% 
% legend
% xlabel("Stimulus")
% ylabel("PDF")
% hold off
% 
% 
% % High confidence report
% subplot(2, 2, 2)
% hold on
% prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
% plot(stimVals, prop_HC, DisplayName="\sigma_m = 0.2", LineWidth=1.5, LineStyle="-")
% legend
% xlabel("Stimulus")
% ylabel("%HC")
% hold off
% 
% 
% stimVals = linspace(-20, 20, 100);
% 
% modelParams.sigma_e                   = 1;
% modelParams.sigma_i                   = 1;
% modelParams.varGain                   = 0.5;
% modelParams.internalNoiseSamplesCnt   = 200;
% modelParams.Cd                        = 0;
% modelParams.Cc                        = 1;
% modelParams.sigma_m                   = 0.2;
% modelParams.sampleCnt                 = 100;
% 
% choicePDFs = getLLhChoice_MCASANDRE(stimVals, modelParams);
% 
% 
% % High confidence report
% subplot(2, 2, 2)
% hold on
% prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
% plot(stimVals, prop_HC, DisplayName="\sigma_m = 0.2", LineWidth=1.5, LineStyle="-")
% legend
% xlabel("Stimulus")
% ylabel("%HC")
% hold off
% 
% 
% stimVals = linspace(-20, 20, 100);
% 
% modelParams.sigma_e                   = 0.1;
% modelParams.sigma_i                   = 1;
% modelParams.varGain                   = 0.5;
% modelParams.internalNoiseSamplesCnt   = 200;
% modelParams.Cd                        = 0;
% modelParams.Cc                        = 1;
% modelParams.sigma_m                   = 0.2;
% modelParams.sampleCnt                 = 100;
% 
% choicePDFs = getLLhChoice_MCASANDRE(stimVals, modelParams);
% 
% 
% % High confidence report
% subplot(2, 2, 2)
% hold on
% prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
% plot(stimVals, prop_HC, DisplayName="\sigma_m = 0.2", LineWidth=1.5, LineStyle="-")
% legend
% xlabel("Stimulus")
% ylabel("%HC")
% hold off
% 
