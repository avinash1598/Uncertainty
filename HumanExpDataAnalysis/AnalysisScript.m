restoredefaultpath
close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')

data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR31.mat'); % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR33.mat'); % Akash

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR32.mat'); % Jiaming
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/ExpScript/CORNFB01.mat');   % Yichao

% formattedData = formatExpData(data);
fltData = data.dat( data.dat.session > 0 , :); %data.dat;
f.dat = fltData; %data.dat; %fltData;
formattedData = formatExpData(f, false, true);
rvOriErr = -90:3:90;

%% Sample mean and dispersion distribution

actualOris      = fltData.stimOri;
actualSpreads   = fltData.stimSpread;
sampleMeans     = fltData.stimSampleMeanOri;
sampleSpreads   = fltData.stimSampleSpread;

figure

% Mean ori distribution
err = sampleMeans - actualOris;
err = mod(err + 90, 180) - 90; % Error between -90 and 90
subplot(2, 2, 1)
histogram(err)
xlabel("Error (deg)")
ylabel("count")
title("Sample mean error dsitribution")

% Sample Spread distribution
err = sampleSpreads - actualSpreads;
subplot(2, 2, 2)
histogram(err)
xlabel("Error (AU)")
ylabel("count")
title("Sample spread error distribution")

x = unique(actualSpreads);
for i=1:numel(x)
    idx = actualSpreads == x(i);

    o1 = actualOris(idx);
    s1 = actualSpreads(idx);
    o2 = sampleMeans(idx);
    s2 = sampleSpreads(idx);
    
    subplot(2, 2, 3)
    hold on
    % Mean ori distribution
    err = o2 - o1;
    err = mod(err + 90, 180) - 90; % Error between -90 and 90
    subplot(2, 2, 3)
    histogram(err, DisplayName="Spread=" + num2str(x(i)))
    xlabel("Error (deg)")
    ylabel("count")
    title("Sample mean error dsitribution")
    legend
    hold off

    % Sample Spread distribution
    err = s2 - s1;
    subplot(2, 2, 4)
    hold on
    histogram(err, DisplayName="Spread=" + num2str(x(i)))
    xlabel("Error (AU)")
    ylabel("count")
    title("Sample spread error distribution")
    legend
    hold off

end

%% Sample stim ori spread by uncertainty level
stimOris             = formattedData.theta_true_all;
sampleStimOri        = formattedData.theta_true_all_S;
conf_reports         = formattedData.confidence_report_all;
n_uncertainty_levels = formattedData.n_uncertainty_levels;

stimOris      = reshape(stimOris, size(stimOris, 1), []);
sampleStimOri = reshape(sampleStimOri, size(sampleStimOri, 1), []);
conf_reports  = reshape(conf_reports, size(conf_reports, 1), []);
oriOffset     = sampleStimOri - stimOris;

figure

subplot(2, 3, 1)
errorbar(1:n_uncertainty_levels, ...
    mean(oriOffset, 2), std(oriOffset, [], 2), 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel("Uncertainty level")
ylabel("Stim - Sample ori")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);

% By HC and LC
oriOffset_HC = oriOffset;
oriOffset_HC(conf_reports == 0) = nan;
oriOffset_LC = oriOffset;
oriOffset_LC(conf_reports == 1) = nan;


subplot(2, 3, 2)
hold on
errorbar(1:n_uncertainty_levels, ...
    mean(oriOffset_HC, 2, 'omitnan'), std(oriOffset, [], 2, 'omitnan'), 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="HC");
errorbar(1:n_uncertainty_levels, ...
    mean(oriOffset_LC, 2, 'omitnan'), std(oriOffset, [], 2, 'omitnan'), 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="LC");

xlabel("Uncertainty level")
ylabel("Stim - Sample ori")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);
legend
hold off

%% Aggregate stats

% Raw err and Confidence (aggregate and by orientation)
figure

orientations               = unique(formattedData.theta_true_all);
resp_err_all               = formattedData.resp_err_all;
confidence_report_all      = formattedData.confidence_report_all;
resp_err_all_flat          = formattedData.resp_err_all(:);
confidence_report_all_flat = formattedData.confidence_report_all(:);

A                          = permute(resp_err_all, [2 1 3:ndims(resp_err_all)]);
resp_err_by_ori            = reshape(A, size(A, 1), []); % is this correct
A                          = permute(confidence_report_all, [2 1 3:ndims(confidence_report_all)]);
confidence_report_by_ori   = reshape(A, size(A, 1), []);

resp_err_all_flat_HC = resp_err_all_flat(confidence_report_all_flat == 1);
resp_err_all_flat_LC = resp_err_all_flat(confidence_report_all_flat == 0);

subplot(3, 3, 1)
histogram(resp_err_all_flat, Normalization="pdf", BinEdges=rvOriErr)
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("All reports, Std: %.2f", std(resp_err_all_flat)))

subplot(3, 3, 2)
r1 = median(resp_err_all_flat_HC) - 5*mad(resp_err_all_flat_HC, 1);
r2 = median(resp_err_all_flat_HC) + 5*mad(resp_err_all_flat_HC, 1);
hold on
histogram(resp_err_all_flat_HC, Normalization="pdf", BinEdges=rvOriErr)

xline(r1, LineStyle="--", LineWidth=1.2);
xline(r2, LineStyle="--", LineWidth=1.2);
hold off
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("HC, Std: %.2f", std(resp_err_all_flat_HC)))

subplot(3, 3, 3)
r1 = median(resp_err_all_flat_LC) - 5*mad(resp_err_all_flat_LC, 1);
r2 = median(resp_err_all_flat_LC) + 5*mad(resp_err_all_flat_LC, 1);
hold on
histogram(resp_err_all_flat_LC, Normalization="pdf", BinEdges=rvOriErr)

xline(r1, LineStyle="--", LineWidth=1.2);
xline(r2, LineStyle="--", LineWidth=1.2);
hold off
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("LC, Std: %.2f", std(resp_err_all_flat_LC)))

propHC = numel(resp_err_all_flat_HC) / numel(resp_err_all_flat);
propLC = numel(resp_err_all_flat_LC) / numel(resp_err_all_flat);

subplot(3, 3, 4)
bar([0, 1], [propHC, propLC], 0.5, DisplayName="data")
xticks([0, 1])
xticklabels({'HC', 'LC'})
ylabel("Proportion")

prop_LC_by_ori = arrayfun(@(i) numel(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 0))/size(resp_err_by_ori,2), 1:size(resp_err_by_ori,1))';
prop_HC_by_ori = arrayfun(@(i) numel(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 1))/size(resp_err_by_ori,2), 1:size(resp_err_by_ori,1))';

subplot(3, 3, 5)
hold on
scatter(orientations, prop_HC_by_ori, DisplayName="HC")
scatter(orientations, prop_LC_by_ori, DisplayName="LC")
plot(orientations, prop_HC_by_ori, DisplayName='HC', HandleVisibility="off")
plot(orientations, prop_LC_by_ori, DisplayName='LC', HandleVisibility="off")
hold off
xlabel("Orientation (deg)")
ylabel("Proportion")
legend

std_err_HC = std(resp_err_all_flat_HC);
std_err_LC = std(resp_err_all_flat_LC);

subplot(3, 3, 6)
bar([0, 1], [std_err_HC, std_err_LC], 0.5, DisplayName="data")
xticks([0, 1])
xticklabels({'HC', 'LC'})
legend
ylabel("std(data)")

% Get LC and HC std per row
std_LC_by_ori = arrayfun(@(i) std(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 0)), 1:size(resp_err_by_ori,1))';
std_HC_by_ori = arrayfun(@(i) std(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 1)), 1:size(resp_err_by_ori,1))';

subplot(3, 3, 7)
hold on
scatter(orientations, std_HC_by_ori, DisplayName='HC')
scatter(orientations, std_LC_by_ori, DisplayName='LC')
plot(orientations, std_HC_by_ori, DisplayName='HC', HandleVisibility="off")
plot(orientations, std_LC_by_ori, DisplayName='LC', HandleVisibility="off")
hold off
xlabel("Oreintation")
ylabel("std(data)")
legend

% Bias 
mean_err   = mean(resp_err_by_ori, 2);
std_val    = std(resp_err_by_ori, [], 2);

subplot(3, 3, 8)
scatter(orientations, mean_err)
hold on
plot(orientations, mean_err)
xlabel("orientation (deg)")
ylabel("Mean error (deg)")
title("Orientation Bias")
hold off

subplot(3, 3, 9)
bar(orientations, std_val);
hold on
xlabel("orientation (deg)")
ylabel("std(data)")
ylim([min(std_val) - 0.1*min(std_val), max(std_val) + 0.1*max(std_val)])
title("Stim dependent std dev")
hold off

%% Response distribution (good if it is uniform)
figure

subplot(2, 3, 1)
histogram(formattedData.theta_resp_all(:), 'BinEdges', 0:5:180);
xlabel("Orientation")
ylabel("Count")
title("Response distribution")

subplot(2, 3, 2)
a = formattedData.theta_resp_all(:);
b = a(formattedData.confidence_report_all(:) == 1);
histogram(b, 'BinEdges', 0:5:180);
xlabel("Orientation")
ylabel("Count")
title("Response distribution (HC)")

subplot(2, 3, 3)
a = formattedData.theta_resp_all(:);
b = a(formattedData.confidence_report_all(:) == 0);
histogram(b, 'BinEdges', 0:5:180);
xlabel("Orientation")
ylabel("Count")
title("Response distribution (LC)")

%% Aggregate stats by uncertainty level

n_uncertainty_levels = formattedData.n_uncertainty_levels;
A                            = resp_err_all;
resp_err_by_level            = reshape(A, size(A, 1), []); % is this correct
A                            = confidence_report_all;
confidence_report_by_level   = reshape(A, size(A, 1), []);

% Raw errors by group - Estimation error (Human subject)
figure
subplot(1, 2, 1)
hold on

for i = 1:n_uncertainty_levels 
    grpOriErr = resp_err_by_level(i, :);
    meanErr = mean(grpOriErr);
    xPts = i + 0.4*(rand(1, numel(grpOriErr)) - 0.5);
    
    scatter(xPts, grpOriErr, 60, 'LineWidth', 0.5);    
    plot([i - 0.35, i + 0.35], [meanErr, meanErr], LineStyle="-", LineWidth=2, Color='black')
end

ylabel("Estimation error (°)", 'FontSize', 14)
xlabel("Uncertainty level", 'FontSize', 14)
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);
yticks([-90, 0, 90])
ylim([-90, 90])

% Improve figure aesthetics
set(gca, 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out', 'Box', 'off')
hold off

% Error by uncertainty
errHC = resp_err_all_flat(confidence_report_all_flat == 1);
errLC = resp_err_all_flat(confidence_report_all_flat == 0);

abErrHC = abs(errHC);
abErrLC = abs(errLC);

means = [mean(abErrLC), mean(abErrHC)];
sems  = [std(abErrLC)/sqrt(numel(abErrLC)), std(abErrHC)/sqrt(numel(abErrHC))]; % No bias in this plot needs to be shown
colors = [0.85 0.33 0.10; 0.47 0.67 0.19]; 

subplot(1, 2, 2)
hold on

labels = ["Uncertain", "Certain"];
for i = 1:2
    b = bar(i+0.5, means(i), 0.6, DisplayName=labels(i)); % 0.6 = bar width
    b.FaceColor = colors(i,:);
    b.EdgeColor = 'none';
end

% Add error bars
errorbar([1.5 2.5], means, sems, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10, HandleVisibility='off');

ylabel("|Estimation Error| (°)", 'FontSize', 14)
yticks([0, 30])
xticks([])
lgd = legend('\color[rgb]{0.85, 0.33, 0.10} Uncertain','\color[rgb]{0.47, 0.67, 0.19} Certain');
lgd.Box = "off";
set(gca, 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out', 'Box', 'off')
ylim([0, 30]) 
hold off

% PDFs by uncertainty
figure
for i=1:n_uncertainty_levels
    
    subplot(2, n_uncertainty_levels/2, i)
    hold on
    
    grpOriErr    = resp_err_by_level(i, :);
    histogram(grpOriErr, rvOriErr, Normalization="pdf");
    
    xlabel("Orientation (deg)")
    ylabel("count")
    title(sprintf("Std %.2f", std(grpOriErr)))
    
    hold off
end

% PDFs by confidence and uncertainty
% Histogram by confidence
figure

confLabels = ["LC", "HC"];
for confVal = [0, 1]
    for i=1:n_uncertainty_levels
        
        subplot(2, n_uncertainty_levels/2, i)
        hold on

        grpOriErr    = resp_err_by_level(i, :);
        reportedConf = confidence_report_by_level(i, :);
        grpOriErr    = grpOriErr(reportedConf == confVal);
        
        histogram(grpOriErr, rvOriErr, Normalization="pdf", DisplayName=confLabels(confVal + 1));
        
        xlabel("Orientation (deg)")
        ylabel("count")
        title(sprintf('C: %.2f D: %.2f T: %.2f ', ...
            formattedData.uncertaintyVals(i, 1), ...
            formattedData.uncertaintyVals(i, 2), ...
            formattedData.uncertaintyVals(i, 3)))
        legend

        hold off
    end
end

% Error and std by confidence
% Bootstrap std
nIter = 5000;
sdHC_boot = zeros(n_uncertainty_levels, nIter);
sdLC_boot = zeros(n_uncertainty_levels, nIter);

for l=1:n_uncertainty_levels
    for i = 1:nIter
        
        HC_errors_ = resp_err_by_level(l, confidence_report_by_level(l, :) == 1);
        LC_errors_ = resp_err_by_level(l, confidence_report_by_level(l, :) == 0);
        
        % fprintf('HC: %d, LC: %d \n', numel(HC_errors_), numel(LC_errors_))

        if numel(HC_errors_) <= numel(LC_errors_) % LC samples greater
            % disp("LC greater")
            % sample from group2 to match size of group1
            x2_sample = randsample(LC_errors_, numel(HC_errors_));
            % fprintf('HC: %d, LC: %d \n', numel(HC_errors_), numel(x2_sample))
            sdHC_boot(l, i) = std(HC_errors_);
            sdLC_boot(l, i) = std(x2_sample);
        else % HC samples greater 
            % disp("HC greater")
            % sample from group2 to match size of group1
            x2_sample = randsample(HC_errors_, numel(LC_errors_));
            % fprintf('HC: %d, LC: %d \n', numel(x2_sample), numel(LC_errors_))
            sdHC_boot(l, i) = std(x2_sample);
            sdLC_boot(l, i) = std(LC_errors_);
        end
    end
end

% x_HC_bootstrapped = mean(sdHC_boot, 2);
% x_LC_bootstrapped = mean(sdLC_boot, 2);

figure

x = mean(resp_err_by_level, 2);
% x_m = median(resp_err_by_level, 2);
y = std(resp_err_by_level, 0, 2);
% y = mad(resp_err_by_level, 1, 2);

HC_idx = confidence_report_by_level == 1;
LC_idx = confidence_report_by_level == 0;

resp_HC = resp_err_by_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_level;
resp_LC(~LC_idx) = NaN;

x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');
% y_HC = mad(resp_HC, 1, 2);

x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');
% y_LC = mad(resp_LC, 1, 2);

subplot(3, 3, 1)
errorbar(1:n_uncertainty_levels, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6);

xlabel("Uncertainty level")
ylabel("Std(data)")

subplot(3, 3, 2)

% Behavioral variability
plot(1:n_uncertainty_levels, y, LineWidth=2);
xlabel("Uncertainty level")
ylabel("Std(data)")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);

subplot(3, 3, 3)

% Behavioral variability
plot(1:n_uncertainty_levels, y_HC,  LineWidth=2, DisplayName="HC");
% plot(y, y_HC,  LineWidth=2, DisplayName="HC");
hold on
plot(1:n_uncertainty_levels, y_LC,  LineWidth=2, DisplayName="LC");
% plot(y, y_LC,  LineWidth=2, DisplayName="LC");
legend
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);
xlabel("Uncertainty level")
ylabel("Std(data)")
hold off

subplot(3, 3, 4)
n_trials_by_uncertainty_HC =  sum(confidence_report_by_level, 2);
n_trials_by_uncertainty_LC =  sum(~confidence_report_by_level, 2);

hold on
bar((1:n_uncertainty_levels) - 0.125, n_trials_by_uncertainty_HC, 0.25, DisplayName="HC")
bar((1:n_uncertainty_levels) + 0.125, n_trials_by_uncertainty_LC, 0.25, DisplayName="LC")
xlabel("Uncertainty level")
ylabel("Num trials")
xticks(1:n_uncertainty_levels)
xticklabels(1:n_uncertainty_levels)
legend
hold off


% % Standard error instead of standard deviation
% y_    = std(resp_err_by_level, 0, 2) ./ sqrt(sum(~isnan(resp_err_by_level), 2));
% y_HC_ = std(resp_HC, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(resp_HC), 2));
% y_LC_ = std(resp_LC, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(resp_LC), 2));

% MAD of standard deviation
y_    = mad(resp_err_by_level, 1, 2);
y_HC_ = mad(resp_HC, 1, 2);
y_LC_ = mad(resp_LC, 1, 2);

subplot(3, 3, 5)

% Behavioral variability
plot(1:n_uncertainty_levels, y_, LineWidth=2);
xlabel("Uncertainty level")
% ylabel("SEM(data)")
ylabel("MAD")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);

% Behavioral variability
subplot(3, 3, 6)
plot(1:n_uncertainty_levels, y_HC_,  LineWidth=2, DisplayName="HC");
hold on
plot(1:n_uncertainty_levels, y_LC_,  LineWidth=2, DisplayName="LC");
xlabel("Uncertainty level")
% ylabel("SEM(data)")
ylabel("MAD")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);
legend
hold off


% Absolute error
x = mean(abs(resp_err_by_level), 2);
y = std(abs(resp_err_by_level), 0, 2);

x_HC = mean(abs(resp_HC), 2, 'omitnan');
y_HC = std(abs(resp_HC), 0, 2, 'omitnan');

x_LC = mean(abs(resp_LC), 2, 'omitnan');
y_LC = std(abs(resp_LC), 0, 2, 'omitnan');

subplot(3, 3, 7)
errorbar(1:n_uncertainty_levels, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6);

xlabel("Uncertainty level")
ylabel("Abs Error (deg)")

subplot(3, 3, 8)
errorbar(1:n_uncertainty_levels, ...
    x_HC, y_HC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="HC");
hold on
errorbar(1:n_uncertainty_levels, ...
    x_LC, y_LC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="LC");

xlabel("Uncertainty level")
ylabel("Abs (error)")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);
legend
hold off

% Behavioral variability - bootstrapped
y_HC_bootstrapped = mean(sdHC_boot, 2);
y_LC_bootstrapped = mean(sdLC_boot, 2);

subplot(3, 3, 9)
plot(1:n_uncertainty_levels, y_HC_bootstrapped,  LineWidth=2, DisplayName="HC");
hold on
plot(1:n_uncertainty_levels, y_LC_bootstrapped,  LineWidth=2, DisplayName="LC");
xlabel("Uncertainty level")
ylabel("Std (bootstrapped data)")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);
legend
hold off


%% Marginal: By  ori
cardinalData = fltData;
cardinalData.cardinalStimOri = abs( floor(fltData.stimOri/90) * 90 - mod(fltData.stimOri, 90) ); %abs(45 - abs(mod(datOK.stimOri, 90) - 45));
tabData = groupsummary(cardinalData, {'cardinalStimOri'}, {'mean', 'std', @(x) mad(x,1), 'numel'}, 'rawOriError');
% [vals, idx] = sort(tabData.std_rawOriError, 'ascend');
[vals, idx] = sort(tabData.fun1_rawOriError, 'ascend');

figure

subplot(3, 4, 1)
x = tabData.cardinalStimOri(idx);
% y = tabData.std_rawOriError(idx);
y = tabData.fun1_rawOriError(idx);

plot(1:numel(x), y, LineWidth=1.5);
xticks(1:numel(x))
xticklabels(x);
xlabel("Orientation")
ylabel("std(raw error)")

tabData = groupsummary(cardinalData, {'cardinalStimOri', 'reportedConf'}, {'mean', 'std',  @(x) mad(x,1), 'numel'}, 'rawOriError');
x_HC = tabData.cardinalStimOri(tabData.reportedConf == 1);
% y_HC = tabData.std_rawOriError(tabData.reportedConf == 1);
y_HC = tabData.fun1_rawOriError(tabData.reportedConf == 1);

x_LC = tabData.cardinalStimOri(tabData.reportedConf == 0);
% y_LC = tabData.std_rawOriError(tabData.reportedConf == 0);
y_LC = tabData.fun1_rawOriError(tabData.reportedConf == 0);

assert(isequal(x_HC, x_LC))

% [~, idx] = ismember(x, x_HC);
x_HC = x_HC(idx);
x_LC = x_LC(idx);
y_HC = y_HC(idx);
y_LC = y_LC(idx);

subplot(3, 4, 2)
hold on
plot(1:numel(x_HC), y_HC, LineWidth=1.5, DisplayName="HC");
plot(1:numel(x_LC), y_LC, LineWidth=1.5, DisplayName="LC");
xticks(1:numel(x_HC))
xticklabels(x_HC);
xlabel("Orientation")
ylabel("std(raw error)")
legend
hold off

subplot(3, 4, 3)
hold on
bar((1:numel(x_HC)) - 0.125, tabData.GroupCount(tabData.reportedConf == 1), 0.25, DisplayName="HC")
bar((1:numel(x_LC)) + 0.125, tabData.GroupCount(tabData.reportedConf == 0), 0.25, DisplayName="LC")
xlabel("Orientation")
ylabel("Trial count")
xticks(1:numel(x_HC))
xticklabels(x_HC);
legend
hold off

% Plot distribution of HC and LC report for each orientation
[G, oriLevels, reportedConf] = findgroups(cardinalData.cardinalStimOri, cardinalData.reportedConf);

errorDist = splitapply(@(x){x}, ...
    cardinalData.rawOriError, G);

tabData = table(oriLevels, reportedConf, errorDist, ...
    'VariableNames', {'cardinalStimOri', 'reportedConf', 'rawOriErrorDist'});

for i = 1:numel(unique(cardinalData.cardinalStimOri))
    ori = tabData.cardinalStimOri(2*i);

    
    tmp = tabData.rawOriErrorDist(2*i - 1);
    errsLC = tmp{1};

    tmp = tabData.rawOriErrorDist(2*i);
    errsHC = tmp{1};
    
    subplot(3, 4, 5 + i - 1)
    hold on
    histogram(errsLC, rvOriErr, DisplayName="LC")
    histogram(errsHC, rvOriErr, DisplayName="HC")
    hold off
    legend
    title(sprintf("Ori: %d", ori))

end

%% Marginal distribution - Alternative method

% Contrast
c_levels                     = formattedData.marginal_contrast_vals;
resp_err_by_c_level          = formattedData.marginal_contrast_resp_err_all;
confidence_report_by_c_level = formattedData.marginal_contrast_confidence_report_all;

HC_idx = confidence_report_by_c_level == 1;
LC_idx = confidence_report_by_c_level == 0;

resp_HC = resp_err_by_c_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_c_level;
resp_LC(~LC_idx) = NaN;

figure
for i=1:numel(c_levels)
    subplot(3, 4, i) 
    hold on
    histogram(resp_err_by_c_level(i, :), rvOriErr)
    xlabel("Error")
    title(sprintf("C: %.2f", c_levels(i)))
    hold off

    subplot(3, 4, i + 2)
    hold on
    histogram(resp_HC(i, :), rvOriErr, DisplayName="HC")
    histogram(resp_LC(i, :), rvOriErr, DisplayName="LC")
    xlabel("Error")
    legend
    title(sprintf("C: %.2f", c_levels(i)))
    hold off

end


% Spread
s_levels                     = formattedData.marginal_spread_vals;
resp_err_by_s_level          = formattedData.marginal_spread_resp_err_all;
confidence_report_by_s_level = formattedData.marginal_spread_confidence_report_all;

HC_idx = confidence_report_by_s_level == 1;
LC_idx = confidence_report_by_s_level == 0;

resp_HC = resp_err_by_s_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_s_level;
resp_LC(~LC_idx) = NaN;

for i=1:numel(s_levels)
    subplot(3, 4, 4 + i) 
    hold on
    histogram(resp_err_by_s_level(i, :), rvOriErr)
    xlabel("Error")
    title(sprintf("S: %.2f", s_levels(i)))
    hold off

    subplot(3, 4, 4 + i + 2)
    hold on
    histogram(resp_HC(i, :), rvOriErr, DisplayName="HC")
    histogram(resp_LC(i, :), rvOriErr, DisplayName="LC")
    xlabel("Error")
    legend
    title(sprintf("S: %.2f", s_levels(i)))
    hold off

end

% Duration
d_levels                     = formattedData.marginal_dur_vals;
resp_err_by_d_level          = formattedData.marginal_dur_resp_err_all;
confidence_report_by_d_level = formattedData.marginal_dur_confidence_report_all;

HC_idx = confidence_report_by_d_level == 1;
LC_idx = confidence_report_by_d_level == 0;

resp_HC = resp_err_by_d_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_d_level;
resp_LC(~LC_idx) = NaN;

for i=1:numel(d_levels)
    subplot(3, 4, 8 + i) 
    hold on
    histogram(resp_err_by_d_level(i, :), rvOriErr)
    xlabel("Error")
    title(sprintf("D: %.2f", d_levels(i)))
    hold off

    subplot(3, 4, 8 + i + 2)
    hold on
    histogram(resp_HC(i, :), rvOriErr, DisplayName="HC")
    histogram(resp_LC(i, :), rvOriErr, DisplayName="LC")
    xlabel("Error")
    legend
    title(sprintf("D: %.2f", d_levels(i)))
    hold off

end


% 1: Contrast =========================================
figure

c_levels                     = formattedData.marginal_contrast_vals;
resp_err_by_c_level          = formattedData.marginal_contrast_resp_err_all;
confidence_report_by_c_level = formattedData.marginal_contrast_confidence_report_all;

x = mean(resp_err_by_c_level, 2);
y = std(resp_err_by_c_level, 0, 2);
y_m = mad(resp_err_by_c_level, 1, 2);

HC_idx = confidence_report_by_c_level == 1;
LC_idx = confidence_report_by_c_level == 0;

resp_HC = resp_err_by_c_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_c_level;
resp_LC(~LC_idx) = NaN;

% x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');
y_HC_m = mad(resp_HC, 1, 2);

% x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');
y_LC_m = mad(resp_LC, 1, 2);

subplot(3, 5, 1)
errorbar(1:numel(c_levels), ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xticks(1:numel(c_levels))
xticklabels(c_levels);
xlabel("Contrast")
ylabel("Mean err")

subplot(3, 5, 2)

% Behavioral variability
plot(1:numel(c_levels), y, LineWidth=2);
xticks(1:numel(c_levels))
xticklabels(c_levels);
xlabel("Contrast")
ylabel("Std(data)")

subplot(3, 5, 3)

% Behavioral variability
plot(1:numel(c_levels), y_HC,  LineWidth=2, DisplayName="HC");
hold on
plot(1:numel(c_levels), y_LC,  LineWidth=2, DisplayName="LC");
legend
xticks(1:numel(c_levels))
xticklabels(c_levels);
xlabel("Contrast")
ylabel("Std(data)")
hold off

subplot(3, 5, 4)

% Behavioral variability
plot(1:numel(c_levels), y_m, LineWidth=2);
xticks(1:numel(c_levels))
xticklabels(c_levels);
xlabel("Contrast")
ylabel("MAD")

subplot(3, 5, 5)

% Behavioral variability
plot(1:numel(c_levels), y_HC_m,  LineWidth=2, DisplayName="HC");
hold on
plot(1:numel(c_levels), y_LC_m,  LineWidth=2, DisplayName="LC");
legend
xticks(1:numel(c_levels))
xticklabels(c_levels);
xlabel("Contrast")
ylabel("MAD")
hold off

% 2: spread =========================================

s_levels                     = formattedData.marginal_spread_vals;
resp_err_by_s_level          = formattedData.marginal_spread_resp_err_all;
confidence_report_by_s_level = formattedData.marginal_spread_confidence_report_all;

x = mean(resp_err_by_s_level, 2);
y = std(resp_err_by_s_level, 0, 2);
y_m = mad(resp_err_by_s_level, 1, 2);

HC_idx = confidence_report_by_s_level == 1;
LC_idx = confidence_report_by_s_level == 0;

resp_HC = resp_err_by_s_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_s_level;
resp_LC(~LC_idx) = NaN;

% x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');
y_HC_m = mad(resp_HC, 1, 2);

% x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');
y_LC_m = mad(resp_LC, 1, 2);

subplot(3, 5, 6)
errorbar(1:numel(s_levels), ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xticks(1:numel(s_levels))
xticklabels(s_levels);
xlabel("Spread")
ylabel("Mean err")

subplot(3, 5, 7)

% Behavioral variability
plot(1:numel(s_levels), y, LineWidth=2);
xticks(1:numel(s_levels))
xticklabels(s_levels);
xlabel("Spread")
ylabel("Std(data)")

subplot(3, 5, 8)

% Behavioral variability
plot(1:numel(s_levels), y_HC,  LineWidth=2, DisplayName="HC");
hold on
plot(1:numel(s_levels), y_LC,  LineWidth=2, DisplayName="LC");
legend
xticks(1:numel(s_levels))
xticklabels(s_levels);
xlabel("Spread")
ylabel("Std(data)")
hold off

subplot(3, 5, 9)

% Behavioral variability
plot(1:numel(s_levels), y_m, LineWidth=2);
xticks(1:numel(s_levels))
xticklabels(s_levels);
xlabel("Spread")
ylabel("MAD")

subplot(3, 5, 10)

% Behavioral variability
plot(1:numel(s_levels), y_HC_m,  LineWidth=2, DisplayName="HC");
hold on
plot(1:numel(s_levels), y_LC_m,  LineWidth=2, DisplayName="LC");
legend
xticks(1:numel(s_levels))
xticklabels(s_levels);
xlabel("Spread")
ylabel("MAD")
hold off

% 3: Duration =========================================

d_levels                     = formattedData.marginal_dur_vals;
resp_err_by_d_level          = formattedData.marginal_dur_resp_err_all;
confidence_report_by_d_level = formattedData.marginal_dur_confidence_report_all;

x = mean(resp_err_by_d_level, 2);
y = std(resp_err_by_d_level, 0, 2);
y_m = mad(resp_err_by_d_level, 1, 2);

HC_idx = confidence_report_by_d_level == 1;
LC_idx = confidence_report_by_d_level == 0;

resp_HC = resp_err_by_d_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_d_level;
resp_LC(~LC_idx) = NaN;

% x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');
y_HC_m = mad(resp_HC, 1, 2);

% x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');
y_LC_m = mad(resp_LC, 1, 2);

subplot(3, 5, 11)

errorbar(1:numel(d_levels), ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xticks(1:numel(d_levels))
xticklabels(d_levels);
xlabel("Duration")
ylabel("Mean err")

subplot(3, 5, 12)

% Behavioral variability
plot(1:numel(d_levels), y, LineWidth=2);
xticks(1:numel(d_levels))
xticklabels(d_levels);
xlabel("Duration")
ylabel("Std(data)")

subplot(3, 5, 13)

% Behavioral variability
plot(1:numel(d_levels), y_HC,  LineWidth=2, DisplayName="HC");
hold on
plot(1:numel(d_levels), y_LC,  LineWidth=2, DisplayName="LC");
legend
xticks(1:numel(d_levels))
xticklabels(d_levels);
xlabel("Duration")
ylabel("Std(data)")
hold off

subplot(3, 5, 14)

% Behavioral variability
plot(1:numel(d_levels), y_m, LineWidth=2);
xticks(1:numel(d_levels))
xticklabels(d_levels);
xlabel("Duration")
ylabel("MAD")

subplot(3, 5, 15)

% Behavioral variability
plot(1:numel(d_levels), y_HC_m,  LineWidth=2, DisplayName="HC");
hold on
plot(1:numel(d_levels), y_LC_m,  LineWidth=2, DisplayName="LC");
legend
xticks(1:numel(d_levels))
xticklabels(d_levels);
xlabel("Duration")
ylabel("MAD")
hold off


%%
figure

confLabels = ["LC", "HC"];

for i = 1:2
    for confVal = [0 1]

        % contrast
        subplot(2, 3, (i-1)*3 + 1)
        hold on
        
        grpOriErr    = resp_err_by_c_level(i, :);
        reportedConf = confidence_report_by_c_level(i, :);
        grpOriErr    = grpOriErr(reportedConf == confVal);
        
        histogram(grpOriErr, rvOriErr, Normalization="pdf", DisplayName=confLabels(confVal + 1));
        
        xlabel("Err (deg)")
        ylabel("count")
        title(sprintf('Contrast: %.2f ', c_levels(i)))
        legend
        
        hold off


        % Spread
        subplot(2, 3, (i-1)*3 + 2)
        hold on
        
        grpOriErr    = resp_err_by_s_level(i, :);
        reportedConf = confidence_report_by_s_level(i, :);
        grpOriErr    = grpOriErr(reportedConf == confVal);
        
        histogram(grpOriErr, rvOriErr, Normalization="pdf", DisplayName=confLabels(confVal + 1));
        
        xlabel("Err (deg)")
        ylabel("count")
        title(sprintf('Spread: %.2f ', s_levels(i)))
        legend
        
        hold off

        % Duration
        subplot(2, 3, (i-1)*3 + 3)
        hold on
        
        grpOriErr    = resp_err_by_d_level(i, :);
        reportedConf = confidence_report_by_d_level(i, :);
        grpOriErr    = grpOriErr(reportedConf == confVal);
        
        histogram(grpOriErr, rvOriErr, Normalization="pdf", DisplayName=confLabels(confVal + 1));
        
        xlabel("Err (deg)")
        ylabel("count")
        title(sprintf('Duration: %.2f ', d_levels(i)))
        legend
        
        hold off
    end
end


%% Bias and orientation dependent std dev
bias = formattedData.bias;
orientations = formattedData.orientations;

figure
plot(orientations, bias)
xlabel("Orientation")
ylabel("Mean error")
title("Orientation bias")


figure

stdByOri = formattedData.stdByOri;
madByOri = formattedData.madByOri;

for i = 1:n_uncertainty_levels
    subplot(1, 2, 1)
    hold on
    plot(orientations, stdByOri(i, :), LineWidth=1.5, DisplayName=""+i)
    xlabel("Orientation")
    ylabel("Std dev")
    legend
    hold off
    
    subplot(1, 2, 2)
    hold on
    plot(orientations, madByOri(i, :), LineWidth=1.5, DisplayName=""+i)
    xlabel("Orientation")
    ylabel("MAD")
    legend
    hold off
end


% %% Stats split by orientation
% 
% figure(9)
% 
% orientations               = unique(formattedData.theta_true_all);
% n_uncertainty_levels       = formattedData.n_uncertainty_levels;
% n_orientations             = numel(orientations);
% resp_err_all               = formattedData.resp_err_all;
% confidence_report_all      = formattedData.confidence_report_all;
% 
% for i=1:n_uncertainty_levels
%     for j=1:n_orientations
% 
%         errs_         = resp_err_all(i, j, :);
%         conf_reports_ = confidence_report_all(i, j, :);
%         
%         countHC_ = sum(conf_reports_ == 1);
%         countLC_ = sum(conf_reports_ == 0);
%         
%         subplot(n_uncertainty_levels, n_orientations, n_orientations*(i - 1) + j)
%         bar([0, 1], [countHC_, countLC_])
%         xticks([0 1])
%         xticklabels(["HC", "LC"]);
%         ylim([0 25])
% 
%         if j == 1
%             ylabel("count")
%         end
% 
%         if i == 1
%             title(orientations(j))
%         end
%     end
% end
% 
% figure(10)
% 
% for i=1:n_uncertainty_levels
%     for j=1:n_orientations
% 
%         errs_         = resp_err_all(i, j, :);
%         conf_reports_ = confidence_report_all(i, j, :);
%         
%         stdHC_ = std( errs_(conf_reports_ == 1) );
%         stdLC_ = std( errs_(conf_reports_ == 0) );
% 
%         % stdHC_ = mean( errs_(conf_reports_ == 1) );
%         % stdLC_ = mean( errs_(conf_reports_ == 0) );
%         
%         subplot(n_uncertainty_levels, n_orientations, n_orientations*(i - 1) + j)
%         bar([0, 1], [stdHC_, stdLC_])
%         xticks([0 1])
%         xticklabels(["HC", "LC"]);
%         % ylim([0 60])
%         
%         if j == 1
%             ylabel("std (data)")
%         end
%         
%         if i == 1
%             title(orientations(j))
%         end
%     end
% end
% 


