
%% fit HC and LC error distribution for each uncertainty level
close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')

modelData = load('modelContOriData_cov.mat');
grpOriErr = modelData.data.err; 
confReport   = modelData.data.confReport;

orientations = unique( modelData.data.stimOri(:) );
oreintations = orientations';

n_theta              = numel(orientations);
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr_reshaped = reshape(grpOriErr, n_uncertainty_levels, []);

grpOriErr = reshape(grpOriErr, n_uncertainty_levels, []);
confReport = reshape(confReport, n_uncertainty_levels, []);
rvOriErr     = -90:0.1:90;

param_b              = rand(1, n_uncertainty_levels); % std(grpOriErr_reshaped, [], 2)';  
param_a              = rand(1, n_uncertainty_levels); % 0.5*param_b;                    
param_biasAmp        = rand;                           
param_scale          = rand;
param_sigma_meta     = rand;
param_Cc             = rand; 

params = [param_b param_a param_biasAmp param_scale param_sigma_meta param_Cc];

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

% TODO: in actual data check when this is NaN.

metaData.orientations = oreintations;
metaData.rvOriErr     = rvOriErr;
metaData.pdf_stim_HC  = pdf_stim_HC;
metaData.pdf_stim_LC  = pdf_stim_LC;

% figure 
% for i=1:n_uncertainty_levels
% 
%     subplot(4, 5, i)
%     y = pdf_stim_LC(i, :);
%     plot(metaData.rvOriErr, y(:), LineWidth=1.5, DisplayName="LC");
%     hold on
%     y = pdf_stim_HC(i, :);
%     plot(metaData.rvOriErr, y(:), LineWidth=1.5, DisplayName="HC");
%     hold off
%     xline(0, LineStyle="--")
%     ylim([0, 1])
%     xlim([-7, 7])
%     xlabel("Error (deg)")
%     ylabel("count")
%     title("All orientations")
% 
% end

%% Fit model
nParams = numel(params); 

% Objective function
objFun = @(x) minimizeError(x, grpOriErr, metaData);

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

opt_param_b              = optimalValues(1:n_uncertainty_levels);
opt_param_a              = optimalValues(n_uncertainty_levels+1:2*n_uncertainty_levels);
opt_param_biasAmp        = optimalValues(2*n_uncertainty_levels + 1);
opt_param_scale          = optimalValues(2*n_uncertainty_levels + 2);
opt_param_sigma_meta     = optimalValues(2*n_uncertainty_levels + 3);
opt_param_Cc             = optimalValues(2*n_uncertainty_levels + 4);

% %% Plot fit result 
% 
% figure 
% for i=1:n_uncertainty_levels
% 
%     modelParams.b                   = opt_param_b(i);
%     modelParams.a                   = opt_param_a(i);
%     modelParams.biasAmp             = opt_param_biasAmp;
%     modelParams.scale               = opt_param_scale;
%     modelParams.Cc                  = opt_param_Cc;
%     modelParams.sigma_meta          = opt_param_sigma_meta;
%     
%     retData = getEstimationsPDF_cov(oreintations, rvOriErr, modelParams);
%    
%     subplot(2, 4, i)
%     hold on
%     
%     y = pdf_stim_LC(i, :);
%     scatter(rvOriErr, y(:), 'filled', DisplayName="LC");
%     plot(rvOriErr, retData.analyticalPDF_LC, HandleVisibility="off", LineWidth=1.5)
% 
%     y = pdf_stim_HC(i, :);
%     scatter(rvOriErr, y(:), 'filled', DisplayName="HC");
%     plot(rvOriErr, retData.analyticalPDF_HC, HandleVisibility="off", LineWidth=1.5)
% 
%     xline(0, LineStyle="--")
%     ylim([0, 1])
%     xlim([-7, 7])
%     xlabel("Error (deg)")
%     ylabel("count")
%     title("All orientations")
% 
%     legend()
%     hold off
% 
% end

%% Loss function for optimization
function loss = minimizeError(params, data, metaData)

nLevels = size(data, 1);

% Params
param_b              = params(1:nLevels);
param_a              = params(nLevels+1:2*nLevels);
param_biasAmp        = params(2*nLevels + 1);
param_scale          = params(2*nLevels + 1);
param_sigma_meta     = params(2*nLevels + 2);
param_Cc             = params(2*nLevels + 3);

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
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    
    retData = getEstimationsPDF_cov(orientations, rvOriErr, modelParams);
    
    currFit_HC(i, :, :) = retData.analyticalPDF_stim_HC;
    currFit_LC(i, :, :) = retData.analyticalPDF_stim_LC;
end

lossHC = ( currFit_HC - targetPDF_HC ).^2;
lossLC = ( currFit_LC - targetPDF_LC ).^2;

loss = sum( lossHC + lossLC , 'all'); % , 'omitnan'
% loss = mean( lossHC + lossLC , 'all', 'omitnan');

% disp(sum(isnan(targetPDF_HC(:))))
% disp(sum(isnan(targetPDF_LC(:))))
% disp(loss)

end

