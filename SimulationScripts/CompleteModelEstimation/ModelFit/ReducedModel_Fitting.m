
%% fit HC and LC error distribution for each uncertainty level
close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')

modelData = load('modelContOriData.mat');
grpOriErr = modelData.data.err; 
confReport   = modelData.data.confReport;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr = reshape(grpOriErr, n_uncertainty_levels, []);
confReport = reshape(confReport, n_uncertainty_levels, []);
rvOriErr     = -90:0.1:90;

param_sigma_s        = std(grpOriErr, [], 2)';  % rand(1, n_uncertainty_levels); % Choose b such that average noise level ranges from low to high (relative to internal noise level)
param_shape          = rand;
param_scale          = rand;
param_sigma_meta     = rand;
param_Cc             = rand; 

params = [param_sigma_s param_shape param_scale param_sigma_meta param_Cc];

% Get PDFs from data for HC and LC
pdf_stim_LC = zeros( n_uncertainty_levels, numel(rvOriErr) );
pdf_stim_HC = zeros( n_uncertainty_levels, numel(rvOriErr) );

for i=1:n_uncertainty_levels

    cR = confReport(i, :);
    dataHC = grpOriErr(i, cR == 1);
    dataLC = grpOriErr(i, cR == 0);
    
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

HC_idx = confReport == 1;
LC_idx = confReport == 0;

resp_HC = grpOriErr;
resp_HC(~HC_idx) = NaN;

resp_LC = grpOriErr;
resp_LC(~LC_idx) = NaN;

std_HC = std(resp_HC, 0, 2, 'omitnan');
std_LC = std(resp_LC, 0, 2, 'omitnan');

metaData.rvOriErr     = rvOriErr;
metaData.pdf_stim_HC  = pdf_stim_HC;
metaData.pdf_stim_LC  = pdf_stim_LC;
metaData.targetStds   = std( grpOriErr, [], 2 )';
metaData.std_HC       = std_HC';
metaData.std_LC       = std_LC';

figure 
for i=1:n_uncertainty_levels

    subplot(4, 5, i)
    y = pdf_stim_LC(i, :);
    plot(metaData.rvOriErr, y(:), LineWidth=1.5, DisplayName="LC");
    hold on
    y = pdf_stim_HC(i, :);
    plot(metaData.rvOriErr, y(:), LineWidth=1.5, DisplayName="HC");
    hold off
    xline(0, LineStyle="--")
    ylim([0, 1])
    xlim([-7, 7])
    xlabel("Error (deg)")
    ylabel("count")
    title("All orientations")

end

%% Fit model
nParams = numel(params); 

% Objective function
objFun = @(x) minimizeError(x, grpOriErr, metaData);

% Bounds (ga requires finite bounds!)
lb = zeros(size(params));     % same as before
ub = []; % example finite upper bounds

%% Genetic algorithm
% % GA options
% options = optimoptions('ga', ...
%     'Display','iter', ...
%     'UseParallel',false, ...     % turn on parallel if you have Parallel Toolbox
%     'PopulationSize',200, ...
%     'MaxGenerations',500);
% 
% % Run GA
% [optimalValues, fval] = ga(objFun, nParams, [], [], [], [], lb, ub, [], options);

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

opt_param_sigma_s        = optimalValues(1:n_uncertainty_levels);
opt_param_shape          = optimalValues(n_uncertainty_levels + 1);
opt_param_scale          = optimalValues(n_uncertainty_levels + 2);
opt_param_sigma_meta     = optimalValues(n_uncertainty_levels + 3);
opt_param_Cc             = optimalValues(n_uncertainty_levels + 4);

%% Goodness of parameter and curve fit

% Ground truth
gt_sigma_s = modelData.data.params.sigma_s_reduced_model;
gt_shape = modelData.data.params.shape;
gt_scale = modelData.data.params.scale;
gt_sigma_meta = modelData.data.params.sigma_meta;
gt_Cc = modelData.data.params.Cc;

loss_param_fit = sum( (opt_param_sigma_s - gt_sigma_s).^2 ) + ...
    (opt_param_shape - gt_shape).^2 + ...
    (opt_param_scale - gt_scale).^2 + ...
    (opt_param_sigma_meta - gt_sigma_meta).^2 + ...
    (opt_param_Cc - gt_Cc).^2;

loss_curve_fit = 0;

% Curve fit
for i=1:n_uncertainty_levels

    modelParams.sigma_s             = gt_sigma_s(i);
    modelParams.shape               = gt_shape;   
    modelParams.scale               = gt_scale;
    modelParams.Cc                  = gt_sigma_meta;
    modelParams.sigma_meta          = gt_Cc;
    
    retData_gt = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    
    modelParams.sigma_s             = opt_param_sigma_s(i);
    modelParams.shape               = opt_param_shape;   
    modelParams.scale               = opt_param_scale;
    modelParams.Cc                  = opt_param_Cc;
    modelParams.sigma_meta          = opt_param_sigma_meta;
    
    retData_fit = getEstimatesPDFs_reduced_model(rvOriErr, modelParams);
    
    loss_curve_fit = loss_curve_fit + ( retData_fit.analyticalPDF_LC - retData_gt.analyticalPDF_LC ).^2;

    loss_curve_fit = loss_curve_fit + ( retData_fit.analyticalPDF_HC - retData_gt.analyticalPDF_HC ).^2;
end

loss_curve_fit = sqrt(mean(loss_curve_fit));
loss_param_fit = sqrt(loss_param_fit);

%% 
% Display parameters
for i =1:n_uncertainty_levels
    fprintf("GT: %.4f, Fit: %.4f \n", gt_sigma_s(i), opt_param_sigma_s(i))
end

fprintf("GT: %.4f, Fit: %.4f \n", gt_shape, opt_param_shape)
fprintf("GT: %.4f, Fit: %.4f \n", gt_scale, opt_param_scale)
fprintf("GT: %.4f, Fit: %.4f \n", gt_sigma_meta, opt_param_sigma_meta)
fprintf("GT: %.4f, Fit: %.4f \n", gt_Cc, opt_param_Cc)

%% Plot fit result 
anlytcl_sigma_m = zeros(1, n_uncertainty_levels);
anlytcl_sigma_m_HC = zeros(1, n_uncertainty_levels);
anlytcl_sigma_m_LC = zeros(1, n_uncertainty_levels);

figure 
for i=1:n_uncertainty_levels

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
    % ylim([0, 1])
    xlim([-7, 7])
    xlabel("Error (deg)")
    ylabel("count")
    title("All orientations")
    
    legend()
    hold off
    
end

%% Plot results

figure

mean_sigma_s = 1:n_uncertainty_levels;

x = mean(grpOriErr, 2);
y = std(grpOriErr, 0, 2);

HC_idx = confReport == 1;
LC_idx = confReport == 0;

resp_HC = grpOriErr;
resp_HC(~HC_idx) = NaN;

resp_LC = grpOriErr;
resp_LC(~LC_idx) = NaN;

x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');

x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');

x1 = resp_HC(1, :); valid_idx = ~isnan(x1); x1 = x1(valid_idx);
x2 = resp_LC(1, :); valid_idx = ~isnan(x2); x2 = x2(valid_idx);
x3 = resp_HC(n_uncertainty_levels, :); valid_idx = ~isnan(x3); x3 = x3(valid_idx);
x4 = resp_LC(n_uncertainty_levels, :); valid_idx = ~isnan(x4); x4 = x4(valid_idx);

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

% loss = sum( lossHC + lossLC , 'all');

loss = sum( lossHC +  lossLC , 'all') + ...
    sum( ( curr_sigma_m - targetStds ).^2 ) + ...
    sum( ( curr_sigma_m_HC - targetStds_HC ).^2 ) + ...
    sum( ( curr_sigma_m_LC - targetStds_LC ).^2 );

% disp(loss)

end

