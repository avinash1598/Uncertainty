restoredefaultpath
% close all
clear all

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
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB01.mat'); % Yichao
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB02.mat');   % Jonathan
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB03.mat'); % David

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR32.mat'); % Jiaming

fltData = data.dat( data.dat.session > 0, :);
f.dat = fltData;
formattedData = formatExpData(f, false, true);
rvOriErr = -90:2:90; % Why does NLL decreases as bin becomes very fine?

orientations               = unique(formattedData.theta_true_all);
resp_err_all               = formattedData.resp_err_all;
confidence_report_all      = formattedData.confidence_report_all;
resp_err_all_flat          = formattedData.resp_err_all(:);
confidence_report_all_flat = formattedData.confidence_report_all(:);

n_uncertainty_levels = formattedData.n_uncertainty_levels;
A = resp_err_all; resp_err_by_level = reshape(A, size(A, 1), []); % is this correct
A = confidence_report_all; confidence_report_by_level = reshape(A, size(A, 1), []);

% Never use raw PDF values as probabilities
% Don’t compare across different bin sizes
% 1. The Mathematical DifferenceThe Negative Log-Likelihood is defined as 
% the sum of the log of your probabilities.
% 
% Option 1 ($-\sum \log(P)$): 
% This treats the data as discrete. By multiplying the density by the bin
% width, you are calculating the probability mass of that bin. This is the
% correct way if you want to compare your model to a categorical 
% distribution.
%
% Option 2 ($-\sum \log(f(x))$): This is the standard NLL for continuous 
% variables. In continuous space, we sum the log of the densities, 
% not the probabilities.

% Summary RecommendationStop multiplying by bin size unless you 
% specifically need a discrete probability. Use the raw density ($pdf(x)$)
% for a continuous NLL.Switch from bins to KDE if your goal is to find the
% "true" underlying shape of the data.Use Cross-Validation to pick the bin 
% size/bandwidth. The "best" NLL is the one that stays low when looking at
% data the model hasn't seen yet.

%%
% PDFs by uncertainty
% figure
% for i=1:n_uncertainty_levels
%     
%     subplot(2, ceil( n_uncertainty_levels/2 ), i)
%     hold on
%     
%     grpOriErr    = resp_err_by_level(i, :);
%     histogram(grpOriErr, rvOriErr, Normalization="pdf");
%     
%     xlabel("Orientation (deg)")
%     ylabel("count")
%     title(sprintf("Std %.2f", std(grpOriErr)))
%     
%     hold off
% end

% PDFs by confidence and uncertainty
% Histogram by confidence
figure

% Method 1: KDE
trlProbs = [];
totalNLL_v1 = 0;
% Method 2: binning and counting
totalNLL = 0;
totalNLL_v2 = 0;
trlProbs2 = [];

confLabels = ["LC", "HC"];
for confVal = [0, 1]
    for i=1:n_uncertainty_levels
        
        subplot(3, n_uncertainty_levels, i + n_uncertainty_levels*confVal)
        hold on

        grpOriErr    = resp_err_by_level(i, :);
        reportedConf = confidence_report_by_level(i, :);
        grpOriErr    = grpOriErr(reportedConf == confVal);
        
        % Probability of given confidence level
        pConfReport = sum(reportedConf == confVal) / numel(reportedConf);
        %disp(pConfReport);
        % pConfReport = 0.5;
        % compute for p(conf) of 0.5 as well - this might increase/decrease the NLL
        % though comapred to the baseline
        
        % Method 1
        [f, xi] = ksdensity(grpOriErr, 'Function', 'pdf', 'Support', [-90 90]);
        f = f ./ trapz(xi, f);
        % diff2 = xi(2:end) - xi(1:end-1);
        %disp(diff2(1))
        % trapz(xi, f)
        
        likelihood = calculateLikelihood(f, xi, grpOriErr);
        de = xi(2) - xi(1);
        likelihood = likelihood.*pConfReport; % *de % Don't multiply by de
        nll_ = -log(likelihood);
        totalNLL_v1 = totalNLL_v1 + sum(nll_);
        trlProbs = [trlProbs likelihood];
%         min(grpOriErr(:))
%         max(grpOriErr(:))
%         
        % Method 2
%         binCenters = rvOriErr;
%         binEdges = [ binCenters' - diff(binCenters) binCenters' + diff(binCenters)];
%         binEdges = unique(binEdges);
%         [pdf_v2, edges] = histcounts(grpOriErr, binEdges, 'Normalization', 'pdf');
%         sum(pdf_v2 <= 0)
%         pdf_v2(pdf_v2 <= 0) = eps;
%         binCenters = edges(1:end-1) + diff(edges)/2;
%         likelihood2 = calculateLikelihood(pdf_v2, binCenters, grpOriErr);
%         likelihood2 = likelihood2.*pConfReport;
%         trlProbs2 = [trlProbs2 likelihood2];
        
        % Bin method
        [pdf_v2, edges] = histcounts(grpOriErr, rvOriErr, 'Normalization', 'pdf');
        pdf_v2(pdf_v2 <= 0) = eps; % empty bin replace with small eps
        binCenters = edges(1:end-1) + diff(edges)/2;
        
        [count, ~] = histcounts(grpOriErr, rvOriErr, 'Normalization', 'count');
        %[count, ~] = histcounts(grpOriErr, binEdges, 'Normalization', 'count');
        de = diff(rvOriErr);
        nll = - count.*log(pdf_v2*pConfReport); % *de(1) + eps
        totalNLL = totalNLL + sum(nll(:));
        
        % finer bin size more penalty due to bin size

        %tmp = pdf_v2*pConfReport;
        %sum(tmp(:) == 0)
        
        histogram(grpOriErr, rvOriErr, Normalization="pdf", DisplayName=confLabels(confVal + 1));
        hold on
        plot(xi, f, DisplayName="PDF", LineWidth=1.5)
        %plot(binCenters, pdf_v2, DisplayName="PDF V2", LineWidth=1.5)
        %hold off
        
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

ll = log(trlProbs); % + eps
nll = -sum(ll(:));

ll2 = log(trlProbs2); % + eps
nll2 = -sum(ll2(:));

fprintf('NLL method 1: %.4f \n', nll);
fprintf('NLL method 1: %.4f \n', totalNLL_v1);
fprintf('NLL method 2: %.4f \n', nll2);
fprintf('NLL method 2: %.4f \n', totalNLL);

%% By ori
nitr = 2;
bestCaseNLLs = zeros(1 ,nitr);

for itr=1:nitr

    totalNLL = 0;
    
    for confVal = [0, 1]
        for i=1:n_uncertainty_levels
            for j=1:numel(orientations)
    
                grpOriErr    = squeeze( resp_err_all(i, j, :) );
                reportedConf = squeeze( confidence_report_all(i, j, :) );
                grpOriErr    = grpOriErr(reportedConf == confVal);
                
                pConfReport = sum(reportedConf == confVal) / numel(reportedConf);
                pStim = 1/numel(orientations);
    
                % P(err, conf, stim) = P(err|conf,stim)P(conf|stim)P(stim)
    
                if numel(grpOriErr) == 0
                    continue
                end
                
                %assert(numel(grpOriErr) > 0);
                
                % Method 1
                [f, xi] = ksdensity(grpOriErr, 'Function', 'pdf', 'Support', [-90 90]);
                %f = f ./ trapz(xi, f);
                likelihood = calculateLikelihood(f, xi, grpOriErr);
                likelihood = likelihood.*pConfReport; % *pStim *de % Don't multiply by de
                nll_ = -log(likelihood);
                totalNLL = totalNLL + sum(nll_);
            end
        end
    end
    
    fprintf('NLL by ori: %.4f \n', totalNLL);

    bestCaseNLLs(itr) = totalNLL;
end

figure
histogram(bestCaseNLLs)
xlabel("NLL")
ylabel("count")

% 70982 - 7541
%% Distribution of values
binsizes = [0.01 0.1 0.5 1 2 3 10 20 30];
totalNLL_v1 = zeros(1, numel(binsizes));
totalNLL_v2 = zeros(1, numel(binsizes));

for bidx=1:numel(binsizes)
    rvOriErr = -90:binsizes(bidx):90;

    confLabels = ["LC", "HC"];
    for confVal = [0, 1]
        for i=1:n_uncertainty_levels

            grpOriErr    = resp_err_by_level(i, :);
            reportedConf = confidence_report_by_level(i, :);
            grpOriErr    = grpOriErr(reportedConf == confVal);

            pConfReport = sum(reportedConf == confVal) / numel(reportedConf);
            % pConfReport = 0.5;

            % Method 1
            [f, xi] = ksdensity(grpOriErr, 'Function', 'pdf', 'Support', [-90 90]);
            % f = f ./ trapz(xi, f);
            likelihood = calculateLikelihood(f, xi, grpOriErr);
            likelihood = likelihood.*pConfReport; % *de % Don't multiply by de
            nll_ = -log(likelihood);
            totalNLL_v1(bidx) = totalNLL_v1(bidx) + sum(nll_);
            
            % Method2
            [pdf, ~] = histcounts(grpOriErr, rvOriErr, 'Normalization', 'pdf');
            pdf(pdf <= 0) = eps; % empty bin replace with small eps
            [count, ~] = histcounts(grpOriErr, rvOriErr, 'Normalization', 'count');
            nll_ = - count.*log(pdf*pConfReport); % *de(1) + eps
            totalNLL_v2(bidx) = totalNLL_v2(bidx) + sum(nll_);

        end
    end
end

figure
histogram(totalNLL_v1, DisplayName="KDE")
hold on
histogram(totalNLL_v2, DisplayName="bin method")
hold off
xlabel("NLL")
ylabel("count")
legend


%% Ddistribution by ori

%%
function likelihood = calculateLikelihood(pdf, x_pts, trlData_errors)
    likelihood = interp1(x_pts, pdf, trlData_errors, 'linear'); % 1e-10
end
