%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/Utils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/OptimizationUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/PlotUtils/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/LLScriptsTrialData/')
% addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/LLScriptsUtils/')

% orientations     = linspace(0, 179, 30); %0:10:180; % 
orientations     = 0:15:175; %linspace(0, 179, 10); % Alert!!!! This has impact on quality of analytical solution
ntrials_per_ori  = 2500; %1000;
b                = linspace(0.1, 1.5, 6); % 1.2 % Choose b such that average noise level ranges from low to high (relative to internal noise level)
a                = 0.67.*b; %0.67   % Does a depend upon b? Yes
% biasAmp          = 0; %10       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% shape            = 0.848;
% scale            = 338.1997;
% sigma_meta       = 41.5023;
% Cc               = 0.109; 
% guessRate        = 0;

biasAmp          = 0.5; %0.5; %2; %0.5; %0.5; % problem at 2       % Does bias depend upon uncertainty level? No. This bias level seems okay.
shape            = 2;
scale            = 0.5; %0.5;
sigma_meta       = 0.2;
Cc               = 0.5; %0.7
guessRate        = 0; %0.08;

b = [6.6680 9.8866 13.4934 21.0658 33.0866 37.9527];
a = 0.1626.*b;
biasAmp          = 1.5791; %0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
shape            = 0.0373;
scale            = 76720.7879; %0.5;
sigma_meta       = 31.1465;
Cc               = 0.0503;   
guessRate        = 0.0350;


% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(b);
% n_uncertainty_levels     = numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s_stim = b' + a'*(abs(sind(( orientations - 90 ) )));
% sigma_s_stim = b' + a'*(abs(sind(( 2*orientations ) )));
bias = biasAmp*sind(2*orientations); 

for l=1:uncertainty_levels
    for i = 1:n_theta
        theta_true = orientations(i);   % True orientation
        trials = ntrials_per_ori;
        
        % Step1: Using estimate of this uncertainty, the subject estimates
        % the orientation
        % Internal estimate (Gaussian noise) - Note: this is not wraped
        % In actual behavioral data this will be wraped
        % theta_est = theta_true + sigma_p * randn(trials, 1);
        % Find doubly stochastic theta
        gain = gamrnd(shape, scale, [trials 1]);
        sigma_si_modulated = gain;
        sigma_m_stim = sqrt(sigma_s_stim(l, i).^2 + sigma_si_modulated);
        mean_m_stim = theta_true + bias(i);
        
        % TODO: take into account bias?
        % Wrap the angle at the plotting stage. Note: warapping should be
        % performed according to the true angle.
        theta_est = mean_m_stim + sigma_m_stim .* randn(trials, 1);
        theta_est = mod(theta_est, 180); 
        
        % Simulat guess rate
        guess_tl_idx = randi([1 trials], floor( trials*guessRate ), 1);
        guessOris = 180*rand(numel(guess_tl_idx), 1);
        theta_est(guess_tl_idx) = guessOris;
        
        assert(numel(sigma_m_stim) == trials);
        
        % Step1: Subject first gets an estimate of its uncertainty
        % Subject’s estimate of their uncertainty (meta-uncertainty)
        mu_log = log(sigma_m_stim.^2 ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2));
        sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m_stim.^2)));
        sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
        
        % Confidence variable
        Vc = 1 ./ sigma_hat;
        
        % Store
        theta_true_all(l, i, :)         = theta_true;
        theta_resp_all(l, i, :)         = theta_est;
        confidence_report_all(l, i, :)  = Vc > Cc;
    end
end

% Plot the performance curves
resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error

resp_err_all_reshaped = reshape(resp_err_all, uncertainty_levels, []);
confidence_report_all_reshaped = reshape(confidence_report_all, uncertainty_levels, []);

% Save model data
data.stimOri                = theta_true_all;
data.reportedOri            = theta_resp_all;
data.resp_err_all           = resp_err_all;
data.confidence_report_all  = confidence_report_all;
data.stdByOri               = squeeze( std(resp_err_all, 0, 3) );
data.madByOri               = squeeze( mad(resp_err_all, 1, 3) );
data.orientations           = orientations';


data.params.sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
data.params.b                     = b;
data.params.a                     = a;
data.params.biasAmp               = biasAmp;
data.params.shape                 = shape;
data.params.scale                 = scale;
data.params.sigma_meta            = sigma_meta;
data.params.Cc                    = Cc;
data.params.guessRate             = guessRate;

%% Get analytical solution
%% Gamma samples count
s_cnts    = [5 10 20 30 50 70 100 200];
bin_width =[0.1 0.5 1 1.5 2 2.2 3 5 10];
errors    = zeros(uncertainty_levels, numel(s_cnts), numel(bin_width));

% Bindwidth threshold Depends upon data
% Samples count of 100 is in general good enough

for i=1:uncertainty_levels
    for j=1:numel(bin_width)
        
        rvOriErr = -90:bin_width(j):90; %
        centers = rvOriErr;
        binWidth = mean(diff(centers));
        edges = [centers - binWidth/2, centers(end) + binWidth/2];
        
        grpOriErr = resp_err_all_reshaped(i, :);
        cR = confidence_report_all_reshaped(i, :);
        dataHC = resp_err_all_reshaped(i, cR == 1);
        dataLC = resp_err_all_reshaped(i, cR == 0);
        
        % GT
        [pdf, ~] = histcounts(grpOriErr, ...
                'Normalization', 'pdf', ...
                'BinEdges', edges);
    
        [pdfHC, ~] = histcounts(dataHC, ...
                'Normalization', 'pdf', ...
                'BinEdges', edges);
        
        [pdfLC, ~] = histcounts(dataLC, ...
                'Normalization', 'pdf', ...
                'BinEdges', edges);
    
        for k=1:numel(s_cnts)
        
            modelParams.b                   = b(i);
            modelParams.a                   = a(i);
            modelParams.shape               = shape;
            modelParams.sigma_s             = data.params.sigma_s_reduced_model(i); % It would fit nicely if it were opt
            modelParams.biasAmp             = biasAmp;
            modelParams.scale               = scale;
            modelParams.Cc                  = Cc;
            modelParams.sigma_meta          = sigma_meta;
            modelParams.guessRate           = guessRate;
            modelParams.internalNoiseSamplesCnt = s_cnts(k);
            
            retData = getEstimatesPDFs(orientations, rvOriErr, modelParams, true);
    
            err = (pdf - retData.analyticalPDF).^2 + ...
                (pdfHC - retData.analyticalPDF_HC).^2 + ...
                (pdfLC - retData.analyticalPDF_LC).^2;
    
            errors(i, k, j) = mean( err(:) );
        end
    end
    
end

figure
for i=1:numel(bin_width)
    subplot(2, 5, i)
    plot(s_cnts, squeeze( errors(:, :, i) ), LineWidth=2)
    xlabel("Gamma sample counts")
    ylabel("Mean Error (GT - Analytical)")
    title(sprintf("Bin width: %.2f", bin_width(i)))
    legend
end

figure
for i=1:numel(s_cnts)
    subplot(2, 5, i)
    plot(bin_width, squeeze( errors(:, i, :) ), LineWidth=2)
    xlabel("Bin width")
    ylabel("Mean Error (GT - Analytical)")
    title(sprintf("Samples cnt: %.2f", s_cnts(i)))
    legend
end

% %% bin width
% errors = zeros(uncertainty_levels, 10);
% bin_width=[0.1 0.5 1 1.5 2 2.5 3 5 10 30 45];
% 
% for i=1:uncertainty_levels
%     for j=1:numel(bin_width)
% 
%         rvOriErr = -90:bin_width(j):90; %
%         centers = rvOriErr;
%         binWidth = mean(diff(centers));
%         edges = [centers - binWidth/2, centers(end) + binWidth/2];
%         
%         grpOriErr = resp_err_all_reshaped(i, :);
%         cR = confidence_report_all_reshaped(i, :);
%         dataHC = resp_err_all_reshaped(i, cR == 1);
%         dataLC = resp_err_all_reshaped(i, cR == 0);
%         
%         % GT
%         [pdf, ~] = histcounts(grpOriErr, ...
%                 'Normalization', 'pdf', ...
%                 'BinEdges', edges);
%     
%         [pdfHC, ~] = histcounts(dataHC, ...
%                 'Normalization', 'pdf', ...
%                 'BinEdges', edges);
%     
%         [pdfLC, ~] = histcounts(dataLC, ...
%                 'Normalization', 'pdf', ...
%                 'BinEdges', edges);
% 
%         modelParams.b                   = b(i);
%         modelParams.a                   = a(i);
%         modelParams.sigma_s             = data.params.sigma_s_reduced_model(i); % It would fit nicely if it were opt
%         modelParams.biasAmp             = biasAmp;
%         modelParams.scale               = scale;
%         modelParams.Cc                  = Cc;
%         modelParams.sigma_meta          = sigma_meta;
%         modelParams.guessRate           = guessRate;
%         modelParams.internalNoiseSamplesCnt = 100;
%         
%         retData = getEstimationsPDF_cov(orientations, rvOriErr, modelParams, true);
%         
%         err = (pdf - retData.analyticalPDF).^2 + ...
%             (pdfHC - retData.analyticalPDF_HC).^2 + ...
%             (pdfLC - retData.analyticalPDF_LC).^2;
%         
%         errors(i, j) = mean( err(:) );
%     end
%     
% end
% 
% figure
% plot(bin_width, errors, LineWidth=2)
% xlabel("Bin width")
% ylabel("Mean Error (GT - Analytical)")
% title("Orientation error bin width")
% legend