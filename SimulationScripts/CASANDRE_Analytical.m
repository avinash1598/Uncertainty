clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils/')

stimVals = linspace(-5, 5, 100);
modelParams.sigma_d    = 1;
modelParams.Cd         = 0;
modelParams.Cc         = 1;
modelParams.sigma_m    = 0.2;
modelParams.sampleCnt  = 100;

choicePDFs = getLLhChoice_CASANDRE(stimVals, modelParams);

figure
subplot(2, 2, 1)

hold on

plot(stimVals, choicePDFs(1, :), DisplayName="D0, HC", LineWidth=1.5, LineStyle="-")
plot(stimVals, choicePDFs(2, :), DisplayName="D0, LC", LineWidth=1.5, LineStyle="-")
plot(stimVals, choicePDFs(3, :), DisplayName="D1, HC", LineWidth=1.5, LineStyle="-")
plot(stimVals, choicePDFs(4, :), DisplayName="D1, LC", LineWidth=1.5, LineStyle="-")

legend
xlabel("Stimulus")
ylabel("PDF")
hold off


% High confidence report
subplot(2, 2, 2)
hold on
prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
plot(stimVals, prop_HC, DisplayName="\sigma_m = 0.2", LineWidth=1.5, LineStyle="-")
legend
xlabel("Stimulus")
ylabel("%HC")
hold off

% Psychometric curve
subplot(2, 2, 3)
hold on
% High confidence
psyFnHC = choicePDFs(3,:) ./ (choicePDFs(1,:) + choicePDFs(3,:));
plot(stimVals, psyFnHC, LineWidth=2, DisplayName="HC")
psyFnLC = choicePDFs(4,:) ./ (choicePDFs(2,:) + choicePDFs(4,:));
plot(stimVals, psyFnLC, LineWidth=2, DisplayName="LC")
xlabel("Orientation")
ylabel("Proportion CCW" + ...
    "" + ...
    "")
title("Psychometric function")
legend()
hold off

% stimVals = linspace(-5, 5, 100);
% modelParams.sigma_d    = 1;
% modelParams.Cd         = 0;
% modelParams.Cc         = 1;
% modelParams.sigma_m    = 1;
% modelParams.sampleCnt  = 100;
% 
% choicePDFs = getLLhChoice_CASANDRE(stimVals, modelParams);
% 
% % High confidence report
% subplot(2, 2, 2)
% hold on
% prop_HC = choicePDFs(1, :) + choicePDFs(3, :);
% plot(stimVals, prop_HC, DisplayName="\sigma_m = 1", LineWidth=1.5, LineStyle="-")
% legend
% xlabel("Stimulus")
% ylabel("%HC")
% hold off


