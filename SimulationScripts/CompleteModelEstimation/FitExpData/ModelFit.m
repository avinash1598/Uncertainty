close all
clear all

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

param_sigma_s        = rand(1, uncertainty_levels); %std(resp_err_all, [], 2)';  % Choose b such that average noise level ranges from low to high (relative to internal noise level)
param_shape          = rand;
param_scale          = rand;
param_sigma_meta     = rand;
param_Cc             = rand; 

params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];

HC_idx = confidence_report_all == 1;
LC_idx = confidence_report_all == 0;

resp_HC = resp_err_all;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all;
resp_LC(~LC_idx) = NaN;

std_HC = std(resp_HC, 0, 2, 'omitnan');
std_LC = std(resp_LC, 0, 2, 'omitnan');

metaData.rvOriErr     = rvOriErr;
metaData.pdf_stim_HC  = pdf_stim_HC;
metaData.pdf_stim_LC  = pdf_stim_LC;
metaData.targetStds   = std( resp_err_all, [], 2 )';
metaData.std_HC       = std_HC';
metaData.std_LC       = std_LC';


% figure 
% for i=1:uncertainty_levels
% 
%     subplot(2, 4, i)
%     y = pdf_stim_LC(i, :);
%     plot(rvOriErr, y(:), LineWidth=1.5, DisplayName="LC");
%     hold on
%     y = pdf_stim_HC(i, :);
%     plot(rvOriErr, y(:), LineWidth=1.5, DisplayName="HC");
%     hold off
%     xline(0, LineStyle="--")
%     % ylim([0, 1])
%      % xlim([-7, 7])
%     xlabel("Error (deg)")
%     ylabel("count")
%     title("All orientations")
% 
% end

%% Fit model
nParams = numel(params); 

% Objective function
objFun = @(x) minimizeError(x, resp_err_all, metaData);

% Bounds (ga requires finite bounds!)
lb = zeros(size(params));     % same as before
ub = []; % example finite upper bounds

%% Fmincon
% Optimization options for fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...          % or 'interior-point', 'trust-region-reflective', etc.
    'MaxIterations', 1000, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6);

% Initial guess (required for fmincon)
x0 = params;   % start in the middle of bounds, for example

% Run fmincon
[optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
    [], [], [], [], lb, ub, [], options);

disp('Optimal parameters:');
disp(optimalValues);
disp('Final objective value:');
disp(fval);

opt_param_sigma_s        = optimalValues(1:uncertainty_levels);
opt_param_shape          = optimalValues(uncertainty_levels + 1);
opt_param_scale          = optimalValues(uncertainty_levels + 2);
opt_param_sigma_meta     = optimalValues(uncertainty_levels + 3);
opt_param_Cc             = optimalValues(uncertainty_levels + 4);


%% Plot fit result 
anlytcl_sigma_m = zeros(1, uncertainty_levels);
anlytcl_sigma_m_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_LC = zeros(1, uncertainty_levels);

figure 
for i=1:uncertainty_levels

    modelParams.sigma_s             = opt_param_sigma_s(i);
    modelParams.shape               = opt_param_shape;   
    modelParams.scale               = opt_param_scale;
    modelParams.Cc                  = opt_param_Cc;
    modelParams.sigma_meta          = opt_param_sigma_meta;
    
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

% mean_sigma_s = opt_param_sigma_s; %1:uncertainty_levels;
mean_sigma_s = 1:uncertainty_levels;

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
errorbar(mean_sigma_s, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("\sigma_s(s)")
ylabel("Error")

subplot(2, 3, 2)

% Behavioral variability
scatter(mean_sigma_s, y, "filled");
hold on
plot(mean_sigma_s, anlytcl_sigma_m, LineWidth=1.5);
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(mean_sigma_s, y_HC, "filled", DisplayName="High confidence");
hold on
plot(mean_sigma_s, anlytcl_sigma_m_HC, LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s, y_LC, "filled", DisplayName="Low confidence");
plot(mean_sigma_s, anlytcl_sigma_m_LC, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
legend
hold off


%% Loss function for optimization
function loss = minimizeError(params, data, metaData)

nLevels = size(data, 1);

% Params
param_sigma_s        = params(1:nLevels);
param_shape          = params(nLevels + 1);
param_scale          = params(nLevels + 2);
param_sigma_meta     = params(nLevels + 3);
param_Cc             = params(nLevels + 4);

% Metadata
rvOriErr             = metaData.rvOriErr;
targetPDF_HC         = metaData.pdf_stim_HC;
targetPDF_LC         = metaData.pdf_stim_LC;
targetStds           = metaData.targetStds;
targetStds_HC        = metaData.std_HC;
targetStds_LC        = metaData.std_LC;

currFit_HC = zeros(nLevels, numel(rvOriErr));
currFit_LC = zeros(nLevels, numel(rvOriErr));
curr_sigma_m = zeros(1, nLevels);
curr_sigma_m_HC = zeros(1, nLevels);
curr_sigma_m_LC = zeros(1, nLevels);

for i=1:nLevels
    
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.shape               = param_shape;   
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    
    retData = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    
    currFit_HC(i, :) = retData.analyticalPDF_HC;
    currFit_LC(i, :) = retData.analyticalPDF_LC;
    curr_sigma_m(i)  = retData.E_sigma_m;
    curr_sigma_m_HC(i) = retData.E_sigma_m_HC;
    curr_sigma_m_LC(i) = retData.E_sigma_m_LC;
end

lossHC = ( currFit_HC - targetPDF_HC ).^2;
lossLC = ( currFit_LC - targetPDF_LC ).^2;

% loss = sum( lossHC +  lossLC , 'all'); % + sum( (curr_sigma_m - targetStds).^2 ); % + sum( ( param_sigma_s - targetStds ).^2 ); % + 0.1*sum( (curr_sigma_m - targetStds).^2 );

loss = sum( lossHC +  lossLC , 'all') + ...
    sum( ( curr_sigma_m - targetStds ).^2 ) + ...
    sum( ( curr_sigma_m_HC - targetStds_HC ).^2 ) + ...
    sum( ( curr_sigma_m_LC - targetStds_LC ).^2 );

% if ~isfinite(loss) || isnan(loss)
%     loss = 1e20;
% end

end

