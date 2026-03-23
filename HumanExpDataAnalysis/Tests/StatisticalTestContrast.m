restoredefaultpath
close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR31.mat'); % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR33.mat'); % Akash
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB01.mat'); % Yichao
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB02.mat');   % Jonathan
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB03.mat'); % David

dataFileNames = ["COR31.mat" "COR33.mat" "CORNFB01.mat" "CORNFB02.mat" "CORNFB03.mat"]; %"CORNFB01.mat" "CORNFB02.mat" "CORNFB03.mat"
fltSessionIdx = [0, 0, 0, 0, 2, 1];

L1_errDiffs     = [];
L2_errDiffs     = [];
deltaErrorDiffs = [];

% Contrast marginal
for fIdx=1:numel(dataFileNames)
    
    fileName = dataFileNames(fIdx);
    data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/' + fileName); % David

    % This needs to be done for each subjects
    uniqContrasts = unique( data.dat.stimContrast );
    uniqSpreads   = unique( data.dat.stimSpread );
    uniqDurations = unique( data.dat.stimDur );

    % Contrast
    for s=1:numel(uniqSpreads)
        for d=1:numel(uniqDurations)
            
            fltData = data.dat( data.dat.session > fltSessionIdx(fIdx), :); % TODO: change as per subject
            fltData = fltData(fltData.stimSpread == uniqSpreads(s) & fltData.stimDur == uniqDurations(d), :);
            f.dat = fltData; %data.dat; %fltData;
            formattedData = formatExpData(f, false, true);
            
            if size(formattedData.uncertaintyVals, 1) == 2
                % data needs to be sorted by uncertainty level
                %assertionChecks(formattedData, 'stimContrast')
                reorderedData = reorderThings(formattedData, 'stimContrast');
                assertionChecks(reorderedData, 'stimContrast')
                [arr1, arr2, arr3] = errorMetrics(reorderedData);
    
                deltaErrorDiffs = [deltaErrorDiffs arr3];
                L1_errDiffs     = [L1_errDiffs arr1];
                L2_errDiffs     = [L2_errDiffs arr2];
                
            end
    
        end
    end
end

[h,p,ci,stats] = ttest(deltaErrorDiffs, 0);
[p,h,stats] = signrank(deltaErrorDiffs, 0);

[h,p,ci,stats] = ttest(deltaErrorDiffs, 0, 'Tail', 'right')
[p,h,stats] = signrank(deltaErrorDiffs, 0, 'Tail', 'right')

figure
subplot(2, 2, 1)
hold on
histogram(deltaErrorDiffs, 'BinEdges',-40:5:30) %
xline(0, 'LineStyle',"--")
xlabel("delta (errDiff l2 - errorDiff l1)")
title(sprintf("Multiplicative/Additive test for stim contrast \n pooled across all subjects"))
ylabel("Count")
ylim([0 20])

subplot(2, 2, 2)
hold on
histogram(L1_errDiffs, DisplayName="Uncert. L1") %'BinEdges',0:5:40, 
histogram(L2_errDiffs, DisplayName="Uncert. L2") % 'BinEdges',0:5:40,
xline(0, 'LineStyle',"--")
xlabel("Error Diffs")
ylabel("Count")
legend

%%
function data = reorderThings(formattedData, key)
    if key == "stimContrast" % decreasing -> increasing uncertainty
        cVals = formattedData.uncertaintyVals(:, 1);
        [~, idx] = sort(cVals,'descend');
    end
    
    if key == "stimSpread" % decreasing -> increasing uncertainty
        sVals = formattedData.uncertaintyVals(:, 2);
        [~, idx] = sort(sVals,'ascend');
    end
    
    if key == "stimDur" % decreasing -> increasing uncertainty
        dVals = formattedData.uncertaintyVals(:, 3);
        [~, idx] = sort(dVals,'descend');
    end

    data.resp_err_all          = formattedData.resp_err_all(idx, :, :);
    data.confidence_report_all = formattedData.confidence_report_all(idx, :, :);
    data.n_uncertainty_levels  = formattedData.n_uncertainty_levels;
    data.orientations          = formattedData.orientations;
    data.uncertaintyVals       = formattedData.uncertaintyVals(idx, :);
end

function assertionChecks(formattedData, key)
    if key == "stimContrast" % decreasing -> increasing uncertainty
        cVals = formattedData.uncertaintyVals(:, 1);
        assert(cVals(2) < cVals(1));
    end
    
    if key == "stimSpread" % decreasing -> increasing uncertainty
        sVals = formattedData.uncertaintyVals(:, 2);
        assert(sVals(2) > sVals(1))
    end
    
    if key == "stimDur" % decreasing -> increasing uncertainty
        dVals = formattedData.uncertaintyVals(:, 3);
        assert(dVals(2) < dVals(1))
    end
end


function [L1_errDiffs, L2_errDiffs, deltaErrorDiffs] = errorMetrics(formattedData)
    % compute is delta errordiff at two uncertainty level
    n_uncertainty_levels = formattedData.n_uncertainty_levels;
    n_orientations       = numel( formattedData.orientations );

    assert(n_uncertainty_levels == 2);

    L1_errDiffs = [];
    L2_errDiffs = [];
    deltaErrorDiffs = [];

    for j=1:n_orientations
        
        % Low uncertainty
        errs_1         = formattedData.resp_err_all(1, j, :);
        conf_reports_1 = formattedData.confidence_report_all(1, j, :);
        
        % High uncertainty
        errs_2         = formattedData.resp_err_all(2, j, :);
        conf_reports_2 = formattedData.confidence_report_all(2, j, :);
        
        errHC_1 = errs_1(conf_reports_1 == 1);
        errLC_1 = errs_1(conf_reports_1 == 0);
        
        errHC_2 = errs_2(conf_reports_2 == 1);
        errLC_2 = errs_2(conf_reports_2 == 0);
        
        minContPerCondition = 4;
        randomSampleCnt = min([numel(errHC_1), numel(errLC_1), numel(errHC_2), numel(errLC_2)]); 

        fprintf('%d, %d %d %d \n', numel(errHC_1), numel(errLC_1), numel(errHC_2), numel(errLC_2))
        fprintf('Random sample count %d \n', randomSampleCnt);
        
        if numel(errHC_1) > minContPerCondition && numel(errLC_1) > minContPerCondition ...
            && numel(errHC_2) > minContPerCondition && numel(errLC_2) > minContPerCondition % 4, even with 5 the test is significant
            
            % For each subsample equal number of values and compute the
            % average error over many permutation
            % Uncertainty level 1 LC
            nitr = 100;
            l1_LC_err = 0;
            for itr=1:nitr
                l1_LC_samples = randsample(errLC_1(:), randomSampleCnt);
                l1_LC_err = l1_LC_err + mean( abs(l1_LC_samples) ) ;
            end
            l1_LC_err = l1_LC_err / nitr;

            % Uncertainty level 1 HC
            nitr = 100;
            l1_HC_err = 0;
            for itr=1:nitr
                l1_HC_samples = randsample(errHC_1(:), randomSampleCnt);
                l1_HC_err = l1_HC_err + mean( abs(l1_HC_samples) ) ;
            end
            l1_HC_err = l1_HC_err / nitr;

            % Uncertainty level 2 LC
            nitr = 100;
            l2_LC_err = 0;
            for itr=1:nitr
                l2_LC_samples = randsample(errLC_2(:), randomSampleCnt);
                l2_LC_err = l2_LC_err + mean( abs(l2_LC_samples) ) ;
            end
            l2_LC_err = l2_LC_err / nitr;

            % Uncertainty level 2 HC
            nitr = 100;
            l2_HC_err = 0;
            for itr=1:nitr
                l2_HC_samples = randsample(errHC_2(:), randomSampleCnt);
                l2_HC_err = l2_HC_err + mean( abs(l2_HC_samples) ) ;
            end
            l2_HC_err = l2_HC_err / nitr;
            
            % TODO: later maybe add absolute as well
            errDiff_l1 = l1_LC_err - l1_HC_err; % Error diff at uncertainty level 1
            errDiff_l2 = l2_LC_err - l2_HC_err; % Error diff at uncertainty level 2

%             errDiff_l1 = abs(errDiff_l1); % Error diff at uncertainty level 1
%             errDiff_l2 = abs(errDiff_l2); 

            deltaErrDiff = errDiff_l2 - errDiff_l1;

            L1_errDiffs = [L1_errDiffs errDiff_l1];
            L2_errDiffs = [L2_errDiffs errDiff_l2];
            deltaErrorDiffs = [deltaErrorDiffs deltaErrDiff];

        else
            %errDiff = nan;
        end 
    end
end
