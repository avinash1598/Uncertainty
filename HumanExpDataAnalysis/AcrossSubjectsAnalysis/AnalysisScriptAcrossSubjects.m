restoredefaultpath
% close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR31.mat'); % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR33.mat'); % Akash
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB01.mat'); % Yichao
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB02.mat');   % Jonathan
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB03.mat'); % David

dataFileNames = ["COR31.mat" "COR33.mat" "CORNFB01.mat" "CORNFB02.mat" "CORNFB03.mat"];
fltSessionIdx = [0, 0, 0, 0, 2, 1];

figure
pltIdx = 1;

% Contrast marginal
for fIdx=1:numel(dataFileNames)

    fileName = dataFileNames(fIdx);
    data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/' + fileName); % David

    % This needs to be done for each subjects
    uniqContrasts = unique( data.dat.stimContrast );
    uniqSpreads   = unique( data.dat.stimSpread );
    uniqDurations = unique( data.dat.stimDur );
    
    % For each contrast, spread and dispersion find the STD and count per
    % condition
    
    % Contrast
    for s=1:numel(uniqSpreads)
        for d=1:numel(uniqDurations)
    
            fltData = data.dat( data.dat.session > fltSessionIdx(fIdx), :);
            fltData = fltData(fltData.stimSpread == uniqSpreads(s) & fltData.stimDur == uniqDurations(d), :);
            f.dat = fltData; %data.dat; %fltData;
            formattedData = formatExpData(f, false, true);
            
            if size(formattedData.uncertaintyVals, 1) == 2
                % Do the calculations only if conditions are two
                [std_RawErr, mad_RawErr, count_RawErr, keyVals, mad_RawErr2]  = getMetrics(fltData, 'stimContrast');
                
                subplot(4, 10, pltIdx)
                hold on
                plot(1:numel(keyVals), mad_RawErr2, LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("STD")
                xlabel("contrast")
                title(sprintf("S=%.2f, D=%.2f", uniqSpreads(s), uniqDurations(s)))

                % STD
                subplot(4, 10, 10+ pltIdx)
                hold on
                plot(1:numel(keyVals), std_RawErr(1, :), DisplayName="LC", LineWidth=1.5);
                plot(1:numel(keyVals), std_RawErr(2, :), DisplayName="HC", LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("STD")
                xlabel("contrast")
                %legend

                % MAD
                subplot(4, 10, 20+pltIdx)
                hold on
                plot(1:numel(keyVals), mad_RawErr(1, :), DisplayName="LC", LineWidth=1.5);
                plot(1:numel(keyVals), mad_RawErr(2, :), DisplayName="HC", LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("MAD")
                xlabel("contrast")
                %legend

                % Count
                subplot(4, 10, 30+pltIdx)
                hold on
                bar((1:numel(keyVals)) - 0.125, count_RawErr(1, :), 0.25, DisplayName="LC");
                bar((1:numel(keyVals)) + 0.125, count_RawErr(2, :), 0.25, DisplayName="HC");
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("Count")
                xlabel("contrast")
                legend

                pltIdx = pltIdx + 1;
            end
        end
    end


end


figure
pltIdx = 1;

% Spread marginal
for fIdx=1:numel(dataFileNames)

    fileName = dataFileNames(fIdx);
    data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/' + fileName); % David
    
    % disp(fIdx)

    % This needs to be done for each subjects
    uniqContrasts = unique( data.dat.stimContrast );
    uniqSpreads   = unique( data.dat.stimSpread );
    uniqDurations = unique( data.dat.stimDur );
    
    % For each contrast, spread and dispersion find the STD and count per
    % condition
    
    % Contrast
    for s=1:numel(uniqContrasts)
        for d=1:numel(uniqDurations)
    
            fltData = data.dat( data.dat.session > fltSessionIdx(fIdx), :);
            fltData = fltData(fltData.stimContrast == uniqContrasts(s) & fltData.stimDur == uniqDurations(d), :);
            if numel(fltData) == 0
                continue
            end
            f.dat = fltData; %data.dat; %fltData;
            formattedData = formatExpData(f, false, true);
            
            if size(formattedData.uncertaintyVals, 1) == 2
                % disp(formattedData.uncertaintyVals)
                % Do the calculations only if conditions are two
                [std_RawErr, mad_RawErr, count_RawErr, keyVals, mad_RawErr2]  = getMetrics(fltData, 'stimSpread');
                
                subplot(4, 15, pltIdx)
                hold on
                plot(1:numel(keyVals), mad_RawErr2, LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("STD")
                xlabel("spread")
                title(sprintf("C=%.2f, D=%.2f", uniqContrasts(s), uniqDurations(s)))

                % STD
                subplot(4, 15, 15 + pltIdx)
                hold on
                plot(1:numel(keyVals), std_RawErr(1, :), DisplayName="LC", LineWidth=1.5);
                plot(1:numel(keyVals), std_RawErr(2, :), DisplayName="HC", LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("STD")
                xlabel("spread")
                
                %legend
                
                % MAD
                subplot(4, 15, 30+pltIdx)
                hold on
                plot(1:numel(keyVals), mad_RawErr(1, :), DisplayName="LC", LineWidth=1.5);
                plot(1:numel(keyVals), mad_RawErr(2, :), DisplayName="HC", LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("MAD")
                xlabel("spread")
                %legend

                % Count
                subplot(4, 15, 45+pltIdx)
                hold on
                bar((1:numel(keyVals)) - 0.125, count_RawErr(1, :), 0.25, DisplayName="LC");
                bar((1:numel(keyVals)) + 0.125, count_RawErr(2, :), 0.25, DisplayName="HC");
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("Count")
                xlabel("spread")
                %legend

                pltIdx = pltIdx + 1;
            end
        end
    end


end

%%
figure
pltIdx = 1;
% Duration marginal
for fIdx=1:numel(dataFileNames)

    fileName = dataFileNames(fIdx);
    data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/' + fileName); % David

    % This needs to be done for each subjects
    uniqContrasts = unique( data.dat.stimContrast );
    uniqSpreads   = unique( data.dat.stimSpread );
    uniqDurations = unique( data.dat.stimDur );
    
    %disp(fileName)
    % For each contrast, spread and dispersion find the STD and count per
    % condition
    
    % Contrast
    for s=1:numel(uniqSpreads)
        for c=1:numel(uniqContrasts)
    
            fltData = data.dat( data.dat.session > fltSessionIdx(fIdx), :);
            fltData = fltData(fltData.stimSpread == uniqSpreads(s) & fltData.stimContrast == uniqContrasts(c), :);
            f.dat = fltData; %data.dat; %fltData;
            formattedData = formatExpData(f, false, true);
            
            if size(formattedData.uncertaintyVals, 1) == 2
                %disp(formattedData.uncertaintyVals)
                % Do the calculations only if conditions are two
                [std_RawErr, mad_RawErr, count_RawErr, keyVals, mad_RawErr2]  = getMetrics(fltData, 'stimDur');
                
                subplot(4, 10, pltIdx)
                hold on
                plot(1:numel(keyVals), mad_RawErr2, LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("STD")
                xlabel("Duration")
                title(sprintf("S=%.2f, C=%.2f", uniqSpreads(s), uniqContrasts(s)))

                % STD
                subplot(4, 10, 10+pltIdx)
                hold on
                plot(1:numel(keyVals), std_RawErr(1, :), DisplayName="LC", LineWidth=1.5);
                plot(1:numel(keyVals), std_RawErr(2, :), DisplayName="HC", LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("STD")
                xlabel("Duration")
                
                % MAD
                subplot(4, 10, 20+pltIdx)
                hold on
                plot(1:numel(keyVals), mad_RawErr(1, :), DisplayName="LC", LineWidth=1.5);
                plot(1:numel(keyVals), mad_RawErr(2, :), DisplayName="HC", LineWidth=1.5);
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("MAD")
                xlabel("Duration")
                
                % Count
                subplot(4, 10, 30+pltIdx)
                hold on
                bar((1:numel(keyVals)) - 0.125, count_RawErr(1, :), 0.25, DisplayName="LC");
                bar((1:numel(keyVals)) + 0.125, count_RawErr(2, :), 0.25, DisplayName="HC");
                xticks(1:numel(keyVals))
                xticklabels(keyVals)
                ylabel("Count")
                xlabel("Duration")
                
                pltIdx = pltIdx + 1;
            end
        end
    end


end


function [std_RawErr, mad_RawErr, count_RawErr, keyVals, mad_RawErr2] = getMetrics(data, key)

% Get count, std, and mad % 'stimContrast'
stimSummary = groupsummary(data, {key, 'reportedConf'}, ...
    {'mean', 'std', @(x) mad(x,1)}, 'rawOriError');

stimSummary2 = groupsummary(data, {key}, ...
    {'mean', 'std', @(x) mad(x,1)}, 'rawOriError');

uniqKeys     = unique(stimSummary.(key));
std_RawErr   = zeros(2, numel(uniqKeys)); % confidence x conditions
mad_RawErr   = zeros(2, numel(uniqKeys));
count_RawErr = zeros(2, numel(uniqKeys));

mad_RawErr2 = zeros(1, numel(uniqKeys));

for k=1:numel(uniqKeys) % unique keys
    for c=1:2 % confidence
        T = stimSummary(stimSummary.(key) == uniqKeys(k) & stimSummary.reportedConf == c-1, :);
        std_RawErr(c, k)   = T.std_rawOriError;
        mad_RawErr(c, k)   = T.fun1_rawOriError;
        count_RawErr(c, k) = T.GroupCount;
    end
    
    T = stimSummary2(stimSummary2.(key) == uniqKeys(k), :);
    %mad_RawErr2(k) = T.fun1_rawOriError;
    mad_RawErr2(k) = T.std_rawOriError;
end

keyVals = uniqKeys;

% TODO: sort here
if key == "stimContrast" % decreasing -> increasing uncertainty
    [~, idx] = sort(uniqKeys,'descend');
end

if key == "stimSpread" % decreasing -> increasing uncertainty
    [~, idx] = sort(uniqKeys,'ascend');
end

if key == "stimDur" % decreasing -> increasing uncertainty
    [~, idx] = sort(uniqKeys,'descend');
end

keyVals      = keyVals(idx);
std_RawErr   = std_RawErr(:, idx);
mad_RawErr   = mad_RawErr(:, idx);
count_RawErr = count_RawErr(:, idx);
mad_RawErr2  = mad_RawErr2(idx);

end

