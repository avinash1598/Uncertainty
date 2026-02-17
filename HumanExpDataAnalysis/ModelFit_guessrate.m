close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Scripts/')

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR31.mat'); % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR32.mat'); % Jiaming
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR33.mat');   % Akash
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/ExpScript/CORNFB01.mat');   % Yichao

% Format data
errBins       = -90:1:90;
fltData       = data.dat( data.dat.session > 0 , :); 
f.dat         = fltData; 
formattedData = formatExpData(f, false, false);  % Keep sortByMAD to false (it is set to false in NLL script)

%% Optimization
result = Optimize(formattedData, errBins, "cov");

%%
n_uncertainty_levels     = formattedData.n_uncertainty_levels;
opt_param_sigma_s        = result.x(1:n_uncertainty_levels);
% opt_param_shape          = result.x(n_uncertainty_levels + 1);
opt_param_scale          = result.x(n_uncertainty_levels + 1);
opt_param_sigma_meta     = result.x(n_uncertainty_levels + 2);
opt_param_Cc             = result.x(n_uncertainty_levels + 3);
opt_param_guessrate      = result.x(n_uncertainty_levels + 4);

% Display parameters
for i =1:n_uncertainty_levels
    fprintf("Fit: %.4f \n", opt_param_sigma_s(i))
end

% fprintf("Shape: %.4f \n", opt_param_shape)
fprintf("Scale: %.4f \n", opt_param_scale)
fprintf("sigma_meta: %.4f \n", opt_param_sigma_meta)
fprintf("Cc: %.4f \n", opt_param_Cc)
fprintf("Guess rate: %.4f \n", opt_param_guessrate)

%% Plot result 
% Cov
paramCovModel = result.x;
plotFitResult_gr(formattedData, paramCovModel, "cov");



% %% Optimization
% result = Optimize(formattedData, errBins, "ind");
% 
% %%
% n_uncertainty_levels     = formattedData.n_uncertainty_levels;
% opt_param_sigma_s        = result.x(1:n_uncertainty_levels);
% opt_param_shape          = result.x(n_uncertainty_levels + 1);
% opt_param_scale          = result.x(n_uncertainty_levels + 2);
% opt_param_sigma_meta     = result.x(n_uncertainty_levels + 3);
% opt_param_Cc             = result.x(n_uncertainty_levels + 4);
% opt_param_guessrate      = result.x(n_uncertainty_levels + 5);
% 
% % Display parameters
% for i =1:n_uncertainty_levels
%     fprintf("Fit: %.4f \n", opt_param_sigma_s(i))
% end
% 
% fprintf("Shape: %.4f \n", opt_param_shape)
% fprintf("Scale: %.4f \n", opt_param_scale)
% fprintf("sigma_meta: %.4f \n", opt_param_sigma_meta)
% fprintf("Cc: %.4f \n", opt_param_Cc)
% fprintf("Guess rate: %.4f \n", opt_param_guessrate)
% 
% %% Plot result 
% % Ind
% paramIndModel = result.x;
% plotFitResult_gr(formattedData, paramIndModel, "ind");
% 
