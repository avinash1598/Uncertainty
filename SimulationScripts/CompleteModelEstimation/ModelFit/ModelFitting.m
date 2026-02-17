
%% fit HC and LC error distribution for each uncertainty level
close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')

modelData = load('modelContOriData.mat');
grpOriErr = modelData.data.err; 

orientations = unique( modelData.data.stimOri(:) );
oreintations = orientations';

n_theta              = numel(orientations);
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr_reshaped = reshape(grpOriErr, n_uncertainty_levels, []);

param_b              = std(grpOriErr_reshaped, [], 2)';  
param_a              = 0.5*param_b;                    
param_biasAmp        = rand;                           
param_shape          = rand;
param_scale          = rand;
param_sigma_meta     = rand;
param_Cc             = rand; 

params = [param_b param_a param_biasAmp param_shape param_scale param_sigma_meta param_Cc];

rvOriErr     = -90:0.1:90;
confReport   = modelData.data.confReport;

% Get PDFs from data for HC and LC
pdf_stim_LC = zeros( n_uncertainty_levels, n_theta, numel(rvOriErr) );
pdf_stim_HC = zeros( n_uncertainty_levels, n_theta, numel(rvOriErr) );

for i=1:n_uncertainty_levels
    for j=1:n_theta
        
        cR = confReport(i, j, :);
        dataHC = grpOriErr(i, j, cR == 1);
        dataLC = grpOriErr(i, j, cR == 0);
        
        centers = rvOriErr;
        binWidth = mean(diff(centers));
        edges = [centers - binWidth/2, centers(end) + binWidth/2];

        [pdfHC, edges] = histcounts(dataHC, ...
            'Normalization', 'pdf', ...
            'BinEdges', edges);

        [pdfLC, edges] = histcounts(dataLC, ...
            'Normalization', 'pdf', ...
            'BinEdges', edges);
        
        pdf_stim_HC(i, j, :) = pdfHC;
        pdf_stim_LC(i, j, :) = pdfLC;
    end
end

metaData.orientations = oreintations;
metaData.rvOriErr     = rvOriErr;
metaData.pdf_stim_HC  = pdf_stim_HC;
metaData.pdf_stim_LC  = pdf_stim_LC;

% figure 
% for i=1:n_theta
% 
%     subplot(4, 5, i)
%     y = pdf_stim_LC(1, i, :);
%     plot(metaData.rvOriErr, y(:), LineWidth=1.5);
%     xline(0, LineStyle="--")
%     ylim([0, 1])
%     xlim([-7, 7])
%     xlabel("Error (deg)")
%     ylabel("count")
%     title(sprintf("HC: Ori %d", orientations(i)))
% 
% end

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

opt_param_b              = optimalValues(1:n_uncertainty_levels);
opt_param_a              = optimalValues(n_uncertainty_levels+1:2*n_uncertainty_levels);
opt_param_biasAmp        = optimalValues(2*n_uncertainty_levels + 1);
opt_param_shape          = optimalValues(2*n_uncertainty_levels + 2);
opt_param_scale          = optimalValues(2*n_uncertainty_levels + 3);
opt_param_sigma_meta     = optimalValues(2*n_uncertainty_levels + 4);
opt_param_Cc             = optimalValues(2*n_uncertainty_levels + 5);

%% Loss function for optimization
function loss = minimizeError(params, data, metaData)

nLevels = size(data, 1);

% Params
param_b              = params(1:nLevels);
param_a              = params(nLevels+1:2*nLevels);
param_biasAmp        = params(2*nLevels + 1);
param_shape          = params(2*nLevels + 2);
param_scale          = params(2*nLevels + 3);
param_sigma_meta     = params(2*nLevels + 4);
param_Cc             = params(2*nLevels + 5);

% Metadata
orientations         = metaData.orientations;
rvOriErr             = metaData.rvOriErr;
targetPDF_HC         = metaData.pdf_stim_HC;
targetPDF_LC         = metaData.pdf_stim_LC;

currFit_HC = zeros(nLevels, numel(orientations), numel(rvOriErr));
currFit_LC = zeros(nLevels, numel(orientations), numel(rvOriErr));

for i=1:nLevels
    
    modelParams.b                   = param_b(i);
    modelParams.a                   = param_a(i);
    modelParams.biasAmp             = param_biasAmp;
    modelParams.shape               = param_shape;   
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    
    retData = getEstimatesPDFs(orientations, rvOriErr, modelParams);
    
    currFit_HC(i, :, :) = retData.analyticalPDF_stim_HC;
    currFit_LC(i, :, :) = retData.analyticalPDF_stim_LC;
end

lossHC = ( currFit_HC - targetPDF_HC ).^2;
lossLC = ( currFit_LC - targetPDF_LC ).^2;

loss = sum( lossHC + lossLC , 'all');

% disp(loss)

end

