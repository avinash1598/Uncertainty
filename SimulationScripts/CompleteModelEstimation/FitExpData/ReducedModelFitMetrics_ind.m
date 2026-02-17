clear all
close all

addpath('data/')

data = load('ExpDataModelFitMetrics_ind.mat'); 

fitMetrics = data.data.T;

%% Plot parameters 
nParams = 12;
columnNames = fitMetrics.Properties.VariableNames;

fitValues = zeros(1, numel(columnNames)-1);

figure

for i = 1:nParams
    x_ = fitMetrics{:, 1};
    % fltIdx = (x_ == trlCnt);

    p_ = fitMetrics{:, i+1};
    % p_ = p_(fltIdx);
    
    max_val = max(p_);
    min_val = min(p_);

    % actual_param_val = groundTruth(i);
    subplot(4, 3, i)
    
    if i == 12
        BinEdges = 0:0.001:1;
    else
        BinEdges = min_val: (max_val - min_val)/1000 : max_val;
    end

    [pdf, edges] = histcounts(p_, ...
        'Normalization', 'pdf', ...
        'BinEdges', BinEdges);

    [val, idx] = max(pdf);
    fitVal = ( edges(idx) + edges(idx + 1) ) / 2;

    fitValues(i) = fitVal;

    hold on
    
    histogram(p_, BinEdges) %BinEdges=0:0.1:4
    xline(fitVal, LineStyle="--", LineWidth=1.5)
    xlabel("Fit values")
    ylabel("count")
    title(columnNames{i+1})
    hold off
end

%%

addpath('LL_scripts/')

data = load('COR15.mat'); % Changed reward function - affects perceptual variance - more relaxed
% data = load('../Data/COR16.mat'); % orientation dependent reward
% data = load('../Data/COR17.mat'); % ori dependent reward - changed reward function c1 = 5
% COR18 - c1=1, c2=0.3, -3 (all HC)
% COR19 - c1=1, c2=0.3, -0.5 (all HC)

stimOri = data.dat.stimOri;
reportedOri = data.dat.reportedOri;

rawError = reportedOri - stimOri;
rawOriError = mod(rawError + 90, 180) - 90;

data.dat.rawOriError = rawOriError;

selectedRawOriError = data.dat.rawOriError;

%% Prepare data structures
% TODO: Arrange the data in following format nlevels, norientations, ntrialPerOri
% Once I have this datastructure I can just use the existing analysis
% script

[G, contrastLevels, spreadLevels, stimDur] = findgroups(data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur);
grpIdxes           = unique(G);
uncertainty_levels = numel(grpIdxes);
nTrials            = numel(rawOriError(G == grpIdxes(1)));

% Arrange groups in increasing uncertainty level
fitStds = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    grpIdx          = grpIdxes(i);
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    % grpRawErr_Flt   = data.dat.rawOriError(G == grpIdx); % rawOriErrorFlt
    
    fltOriErr = grpRawErr_Flt;
    pd = fitdist(fltOriErr(~isnan(fltOriErr)), 'Normal');
    
    mu = pd.mu;
    sigma = pd.sigma;
    fitStds(i) = sigma;
end

[B_, idx] = sort(fitStds);
sidx = grpIdxes(idx);

% Arrange in structured format - same as the one used for analysis
theta_true_all        = zeros(uncertainty_levels, nTrials);
theta_resp_all        = zeros(uncertainty_levels, nTrials); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, nTrials);
resp_err_all          = zeros(uncertainty_levels, nTrials);

for i=1:uncertainty_levels
    grpIdx          = sidx(i);
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    grpStimOri      = data.dat.stimOri(G == grpIdx);
    grpOriResp      = data.dat.reportedOri(G == grpIdx);
    grpReportedConf = data.dat.reportedConf(G == grpIdx);
    
    theta_true_all(i, :)          = grpStimOri;
    theta_resp_all(i, :)          = grpOriResp;
    confidence_report_all(i, :)   = grpReportedConf;
    resp_err_all(i, :)            = grpRawErr_Flt;
end

%% Get PDF for HC and LC
rvOriErr     = -90:3:90;

% Get PDFs from data for HC and LC
pdf_stim_LC = zeros( uncertainty_levels, numel(rvOriErr) );
pdf_stim_HC = zeros( uncertainty_levels, numel(rvOriErr) );

for i=1:uncertainty_levels

    cR = confidence_report_all(i, :);
    dataHC = resp_err_all(i, cR == 1);
    dataLC = resp_err_all(i, cR == 0);
    
    centers = rvOriErr;
    binWidth = mean(diff(centers));
    edges = [centers - binWidth/2, centers(end) + binWidth/2];

    [pdfHC, edges] = histcounts(dataHC, ...
        'Normalization', 'pdf', ...
        'BinEdges', edges);

    [pdfLC, edges] = histcounts(dataLC, ...
        'Normalization', 'pdf', ...
        'BinEdges', edges);
    
    pdf_stim_HC(i, :) = pdfHC;
    pdf_stim_LC(i, :) = pdfLC;
    
end


%% Plot fit result 
uncertainty_levels = 8;
anlytcl_sigma_m = zeros(1, uncertainty_levels);
anlytcl_sigma_m_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_LC = zeros(1, uncertainty_levels);

figure 
for i=1:uncertainty_levels

    modelParams.sigma_s             = fitValues(i);
    modelParams.shape               = fitValues(uncertainty_levels + 1);
    modelParams.scale               = fitValues(uncertainty_levels + 2);
    modelParams.sigma_meta          = fitValues(uncertainty_levels + 3);
    modelParams.Cc                  = fitValues(uncertainty_levels + 4);
    
    retData = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    
    anlytcl_sigma_m(i)    = retData.E_sigma_m;
    anlytcl_sigma_m_HC(i) = retData.E_sigma_m_HC;
    anlytcl_sigma_m_LC(i) = retData.E_sigma_m_LC;
    
    subplot(2, 4, i)
    hold on
    
    y = pdf_stim_LC(i, :);
    scatter(rvOriErr, y(:), 'filled', DisplayName="LC");
    plot(rvOriErr, retData.analyticalPDF_LC, HandleVisibility="off", LineWidth=1.5)
    
    y = pdf_stim_HC(i, :);
    scatter(rvOriErr, y(:), 'filled', DisplayName="HC");
    plot(rvOriErr, retData.analyticalPDF_HC, HandleVisibility="off", LineWidth=1.5)
    
    xline(0, LineStyle="--")
    %ylim([0, 1])
    %xlim([-7, 7])
    xlabel("Error (deg)")
    ylabel("count")
    title("All orientations")

    legend()
    hold off

end

%% Plot results

figure

% % mean_sigma_s = fitValues(1:uncertainty_levels); % hopefully is fitting is right is is also in increasing uncertainty level
% mean_sigma_s = 1:uncertainty_levels; % increasing uncertainty level
% [val, idx] = sort(mean_sigma_s);

mean_sigma_s = fitValues(1:uncertainty_levels); % hopefully is fitting is right is is also in increasing uncertainty level
idx = 1:uncertainty_levels;

x = mean(resp_err_all, 2);
y = std(resp_err_all, 0, 2);

HC_idx = confidence_report_all == 1;
LC_idx = confidence_report_all == 0;

resp_HC = resp_err_all;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all;
resp_LC(~LC_idx) = NaN;

x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');

x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');

x1 = resp_HC(1, :); valid_idx = ~isnan(x1); x1 = x1(valid_idx);
x2 = resp_LC(1, :); valid_idx = ~isnan(x2); x2 = x2(valid_idx);
x3 = resp_HC(uncertainty_levels, :); valid_idx = ~isnan(x3); x3 = x3(valid_idx);
x4 = resp_LC(uncertainty_levels, :); valid_idx = ~isnan(x4); x4 = x4(valid_idx);

subplot(2, 3, 1)
errorbar(mean_sigma_s(idx), ...
    x(idx), y(idx), 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

% xlabel("\sigma_s(s)")
xlabel("\sigma_s")
ylabel("Error")
title(sprintf("points arrnaged in order of \nincreasing uncertainty"))

subplot(2, 3, 2)

% Behavioral variability
scatter(mean_sigma_s(idx), y(idx), "filled");
hold on
plot(mean_sigma_s(idx), anlytcl_sigma_m(idx), LineWidth=1.5);
% xlabel("\sigma_s(s) (sensory noise)")
xlabel("\sigma_s (sensory noise)")
ylabel("\sigma_m(s) (measurement noise)")
title(sprintf("points arrnaged in order of \nincreasing uncertainty"))

hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(mean_sigma_s(idx), y_HC(idx), "filled", DisplayName="High confidence");
hold on
plot(mean_sigma_s(idx), anlytcl_sigma_m_HC(idx), LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s(idx), y_LC(idx), "filled", DisplayName="Low confidence");
plot(mean_sigma_s(idx), anlytcl_sigma_m_LC(idx), LineWidth=1.5, HandleVisibility="off");
% xlabel("\sigma_s(s) (sensory noise)")
% xlabel("increasing sigma_m")
xlabel("\sigma_s (sensory noise)")
ylabel("\sigma_m(s) (measurement noise)")
legend
hold off

% sigma_s should increase with increasing uncertainty
subplot(2, 3, 4)
hold on
scatter(idx, y, DisplayName = "std(data) - \sigma_m")
scatter(idx, mean_sigma_s, DisplayName = "\sigma_s")
xlabel("increasing uncertainty (from std(data))")
legend
hold off


