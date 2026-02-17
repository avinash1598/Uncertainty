%%%
% Questions to think about:
% 1. Does 'a' depend upon 'b'
%
%%%

clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

% sigma_s          = [12.1392   13.9612   19.8033   19.6392   22.2473   26.4453   23.1557   26.5658];
% scale            = 358;
% sigma_meta       = 5.2349;
% Cc               = 0.06; 
% [12.1392   13.9612   19.8033   19.6392   22.2473   26.4453   23.1557   26.5658];

orientations     = 0:15:179; %0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 36; %12; %12; %25;
b                = linspace(1, 2.2, 6); % Note: different minimum noise level (0.1). Choose b such that average noise level ranges from low to high (relative to internal noise level)
a                = 0.67.*b;   % Does a depend upon b? Yes
biasAmp          = 0.5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
scale            = 0.5; %0.5;
sigma_meta       = 0.6;
Cc               = 0.5; 

% orientations     = 0:15:179; %0:10:180; % linspace(0, 180, 18);
% ntrials_per_ori  = 18; %12; %12; %25;
% b                = [12.1392   13.9612   19.8033  22.2473   23.1557   26.5658]; % Note: different minimum noise level (0.1). Choose b such that average noise level ranges from low to high (relative to internal noise level)
% a                = 0*b;   % Does a depend upon b? Yes
% biasAmp          = 5;       % Does bias depend upon uncertainty level? No. This bias level seems okay.
% scale            = 358; %0.5;
% sigma_meta       = 5.2349;
% Cc               = 0.06;

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(b);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

% Simulation loop
% Stimulus dependent sensory noise
sigma_s_stim = b' + a'*(abs(sind(2*orientations)));
bias = biasAmp*sind(2*orientations); 

% sigma_s_stim_conf = b(8) + a(8)*(abs(sind(2*orientations)));
sigma_s_stim_conf = 2.2 + 2.2*1*(abs(sind(2*orientations)));

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
%         shape = sigma_s_stim(l, i);
%         gain = gamrnd(shape, scale, [trials 1]);
%         sigma_si_modulated = gain;
%         sigma_m_stim = sqrt(sigma_s_stim(l, i).^2 + sigma_si_modulated.^2); % For now maybe keep multiplicative noise additive as data seem to show similar trend
%         mean_m_stim = theta_true + bias(i);
        
        % Multiplicative
        shape = sigma_s_stim(l, i).^2 / scale; % divide by scale so that mean is sigma_s
        gain = gamrnd(shape, scale, [trials 1]);
        sigma_m_stim = sqrt( gain );
        % sigma_m_stim = sigma_s_stim(l, i) + zeros(trials, 1);
        mean_m_stim = theta_true + bias(i);
        
        % Poisson noise
        % lambda = sigma_s_stim(l, i);
        % sigma_m_stim = poissrnd(lambda, [trials 1]);
        % mean_m_stim = theta_true + bias(i);
        
        % TODO: take into account bias?
        % Wrap the angle at the plotting stage. Note: warapping should be
        % performed according to the true angle.
        theta_est = mean_m_stim + sigma_m_stim .* randn(trials, 1);
        % Since this is orientation, wrap the angle between 0 and 180
        theta_est = mod(theta_est, 180); 
        % Are the actual distribution near 0 and 180 captured by this in simulation
    
        assert(numel(sigma_m_stim) == trials);
        
        % Step1: Subject first gets an estimate of its uncertainty
        % Subject’s estimate of their uncertainty (meta-uncertainty)
%         mu_log = log(sigma_m_stim.^2 ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2)); % for confidence they just use orientation baseline variance
%         sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m_stim.^2)));
%         sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
        
        % sigma_s_stim_conf

        sigma_m_stim_conf_ = sigma_s_stim_conf(i) + zeros(size(sigma_m_stim));
        mu_log = log(sigma_m_stim_conf_.^2 ./ sqrt(sigma_meta.^2 + sigma_m_stim_conf_.^2)); % for confidence they just use orientation baseline variance
        sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m_stim_conf_.^2)));
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
data.stimOri     = theta_true_all;
data.reportedOri = theta_resp_all;
data.err         = resp_err_all;
data.confReport  = confidence_report_all;

data.params.sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
data.params.b                     = b;
data.params.a                     = a;
data.params.biasAmp               = biasAmp;
data.params.scale                 = scale;
data.params.sigma_meta            = sigma_meta;
data.params.Cc                    = Cc;

% save('modelContOriData_cov.mat', "data")


%% Get analytical solution
anlytcl_sigma_m_stim = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_HC = zeros(1, uncertainty_levels);
anlytcl_sigma_m_stim_LC = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    rvOriErr = -90:0.1:90;
    
    modelParams.b                   = b(i);
    modelParams.a                   = a(i);
    modelParams.biasAmp             = biasAmp;
    modelParams.scale               = scale;
    modelParams.Cc                  = Cc;
    modelParams.sigma_meta          = sigma_meta;
    
    retData = getEstimationsPDF_cov(orientations, rvOriErr, modelParams);
    
    anlytcl_sigma_m_stim(i)    = retData.E_sigma_m;
    anlytcl_sigma_m_stim_HC(i) = retData.E_sigma_m_HC;
    anlytcl_sigma_m_stim_LC(i) = retData.E_sigma_m_LC;
end

%% Response distribution (good if it is uniform)
figure(13)
subplot(2, 3, 1)
histogram(theta_resp_all(:), 'BinEdges', 0:5:180);
xlabel("Orientation")
ylabel("Count")
title("Response distribution")

subplot(2, 3, 2)
a = theta_resp_all(:);
b = a(confidence_report_all(:) == 1);
histogram(b, 'BinEdges', 0:5:180);
xlabel("Orientation")
ylabel("Count")
title("Response distribution (HC)")

subplot(2, 3, 3)
a = theta_resp_all(:);
b = a(confidence_report_all(:) == 0);
histogram(b, 'BinEdges', 0:5:180);
xlabel("Orientation")
ylabel("Count")
title("Response distribution (LC)")

%% Plot results

figure

mean_sigma_s_stim = mean(sigma_s_stim, 2);

x = mean(resp_err_all_reshaped, 2);
y = std(resp_err_all_reshaped, 0, 2);

HC_idx = confidence_report_all_reshaped == 1;
LC_idx = confidence_report_all_reshaped == 0;

resp_HC = resp_err_all_reshaped;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all_reshaped;
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
errorbar(mean_sigma_s_stim, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("\sigma_s(s)")
ylabel("Error")

subplot(2, 3, 2)

% Behavioral variability
scatter(mean_sigma_s_stim, y, "filled");
hold on
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim, LineWidth=1.5);
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(mean_sigma_s_stim, y_HC, "filled", DisplayName="High confidence");
hold on
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_HC, LineWidth=1.5, HandleVisibility="off");
scatter(mean_sigma_s_stim, y_LC, "filled", DisplayName="Low confidence");
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_LC, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
legend
hold off


%% Aggregate stats

% Raw err and Confidence (aggregate and by orientation)
figure

rvOriErr = -10:1:10;
resp_err_all_flat          = resp_err_all(:);
confidence_report_all_flat = confidence_report_all(:);

A                          = permute(resp_err_all, [2 1 3:ndims(resp_err_all)]);
resp_err_by_ori            = reshape(A, size(A, 1), []); % is this correct
A                          = permute(confidence_report_all, [2 1 3:ndims(confidence_report_all)]);
confidence_report_by_ori   = reshape(A, size(A, 1), []);

resp_err_all_flat_HC = resp_err_all_flat(confidence_report_all_flat == 1);
resp_err_all_flat_LC = resp_err_all_flat(confidence_report_all_flat == 0);

subplot(3, 3, 1)
histogram(resp_err_all_flat, Normalization="pdf", BinEdges=rvOriErr)
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("All reports, Std: %.2f", std(resp_err_all_flat)))

subplot(3, 3, 2)
histogram(resp_err_all_flat_HC, Normalization="pdf", BinEdges=rvOriErr)
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("HC, Std: %.2f", std(resp_err_all_flat_HC)))

subplot(3, 3, 3)
histogram(resp_err_all_flat_LC, Normalization="pdf", BinEdges=rvOriErr)
xlabel("Perceptual err (deg)")
ylabel("Count")
title(sprintf("LC, Std: %.2f", std(resp_err_all_flat_LC)))

propHC = numel(resp_err_all_flat_HC) / numel(resp_err_all_flat);
propLC = numel(resp_err_all_flat_LC) / numel(resp_err_all_flat);

subplot(3, 3, 4)
bar([0, 1], [propHC, propLC], 0.5, DisplayName="data")
xticks([0, 1])
xticklabels({'HC', 'LC'})
ylabel("Proportion")

prop_LC_by_ori = arrayfun(@(i) numel(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 0))/size(resp_err_by_ori,2), 1:size(resp_err_by_ori,1))';
prop_HC_by_ori = arrayfun(@(i) numel(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 1))/size(resp_err_by_ori,2), 1:size(resp_err_by_ori,1))';

subplot(3, 3, 5)
hold on
scatter(orientations, prop_HC_by_ori, DisplayName="HC")
scatter(orientations, prop_LC_by_ori, DisplayName="LC")
plot(orientations, prop_HC_by_ori, DisplayName='HC', HandleVisibility="off")
plot(orientations, prop_LC_by_ori, DisplayName='LC', HandleVisibility="off")
hold off
xlabel("Orientation (deg)")
ylabel("Proportion")
legend

std_err_HC = std(resp_err_all_flat_HC);
std_err_LC = std(resp_err_all_flat_LC);

subplot(3, 3, 6)
bar([0, 1], [std_err_HC, std_err_LC], 0.5, DisplayName="data")
xticks([0, 1])
xticklabels({'HC', 'LC'})
legend
ylabel("std(data)")

% Get LC and HC std per row
std_LC_by_ori = arrayfun(@(i) std(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 0)), 1:size(resp_err_by_ori,1))';
std_HC_by_ori = arrayfun(@(i) std(resp_err_by_ori(i, confidence_report_by_ori(i,:) == 1)), 1:size(resp_err_by_ori,1))';

subplot(3, 3, 7)
hold on
scatter(orientations, std_HC_by_ori, DisplayName='HC')
scatter(orientations, std_LC_by_ori, DisplayName='LC')
plot(orientations, std_HC_by_ori, DisplayName='HC', HandleVisibility="off")
plot(orientations, std_LC_by_ori, DisplayName='LC', HandleVisibility="off")
hold off
xlabel("Oreintation")
ylabel("std(data)")
legend

% Bias 
mean_err   = mean(resp_err_by_ori, 2);
std_val    = std(resp_err_by_ori, [], 2);

subplot(3, 3, 8)
scatter(orientations, mean_err)
hold on
plot(orientations, mean_err)
xlabel("orientation (deg)")
ylabel("Mean error (deg)")
title("Orientation Bias")
hold off

subplot(3, 3, 9)
bar(orientations, std_val);
hold on
xlabel("orientation (deg)")
ylabel("std(data)")
ylim([min(std_val) - 0.1*min(std_val), max(std_val) + 0.1*max(std_val)])
title("Stim dependent std dev")
hold off



%% Statistical test
% Is confidence report explained by orientation alone or uncertainty?
x = 1:uncertainty_levels; uncertainty = repmat(x', [1 n_theta ntrials_per_ori]);
certainty_     = confidence_report_all(:);
orientations_  = categorical(theta_true_all(:));
uncertainty_   = categorical(uncertainty(:));
T = table(certainty_, orientations_, uncertainty_, 'VariableNames', {'certainty','orientation', 'uncertaintyLevel'});

Model1 = fitglm(T, ...
    'certainty ~ orientation', ...
    'Distribution', 'binomial', ...
    'Link', 'logit');

% disp(Model1)

Model2 = fitglm(T, ...
    'certainty ~ orientation + uncertaintyLevel', ...
    'Distribution', 'binomial', ...
    'Link', 'logit');
% disp(Model2)

% Extract log-likelihoods
LL1 = Model1.LogLikelihood;
LL2 = Model2.LogLikelihood;

% LR statistic
LR = 2 * (LL2 - LL1);

% Degrees of freedom = # of extra parameters in Model 2
df = Model2.NumCoefficients - Model1.NumCoefficients;

% p-value
p = 1 - chi2cdf(LR, df);

fprintf('\nLikelihood Ratio Test:\n');
fprintf('LR = %.3f, df = %d, p = %.5f\n', LR, df, p);

% p < 0.05 → uncertainty contributes to confidence significantly
% p > 0.05 → confidence explained by orientation alone (orientation-only hypothesis plausible)

%% Are errors split by confidence explained by uncertainty alone of confidence as well?

tmp             = permute(resp_err_all, [2 1 3]);    % original level x ori x trials - converted ori x level x trials
resp_err_by_ori = reshape(tmp, size(tmp, 1), []);    % now 12 × 72
meanErr_by_ori  = mean(resp_err_by_ori, 2);          % Do test with median as well

meanErr_by_ori = reshape(meanErr_by_ori, [1 size(meanErr_by_ori, 1) 1]);
absErr         = (resp_err_all - meanErr_by_ori).^2; %abs( resp_err_all - meanErr_by_ori );
% rawErr = resp_err_all - meanErr_by_ori;

T.abs_error        = absErr(:);
T.log_error        = log(T.abs_error + 1e-6);  % avoid -inf % Approximately normalish
T.orientation      = categorical(T.orientation);
T.uncertaintyLevel = categorical(T.uncertaintyLevel);
T.certainty        = categorical(T.certainty); % HC/LC as categorical

% % M1 = fitlm(T, 'log_error ~ orientation + certainty');
% 
% % disp(M1)
% 
% % M2 = fitlm(T, ...
% %     'log_error ~ orientation + certainty * uncertaintyLevel');
% 
% glme1 = fitglme(T, ...
%     'abs_error ~ certainty + uncertaintyLevel + orientation + certainty * uncertaintyLevel', ... %  + (1|orientation)
%     'Distribution', 'Gamma', 'Link', 'log');
% 
% % disp(glme1)
% 
% % glme2 = fitglme(T, ...
% %     'abs_error ~ certainty + orientation + certainty * orientation', ... %  + (1|orientation)
% %     'Distribution', 'Gamma', 'Link', 'log');
% % 
% % 
% % compare(glme2, glme1)
% % disp(M2)

%%

figure

data = sort(abs(T.abs_error));
n = numel(data);

p = ((1:n)' - 0.5) / n;   % quantile positions

phat = gamfit(data);
theoretical = gaminv(p, phat(1), phat(2));

plot(theoretical, data, 'o'); hold on;
plot(theoretical, theoretical, 'r--');   % identity line
xlabel('Theoretical Gamma Quantiles');
ylabel('Empirical Quantiles');


% % Your data: use absolute error or whatever positive-valued variable
% data = abs(T.abs_error);
% data = data(data > 0);    % remove zeros if present
% 
% % Fit log-normal parameters: mu (mean of log), sigma (std of log)
% phat = lognfit(data);   % phat = [mu, sigma]
% 
% % Sort empirical data
% data_sorted = sort(data);
% n = numel(data_sorted);
% 
% % Empirical quantile positions
% p = ((1:n)' - 0.5) / n;
% 
% % Compute theoretical log-normal quantiles
% theoretical = logninv(p, phat(1), phat(2));
% 
% % Plot Q–Q
% figure; 
% plot(theoretical, data_sorted, 'o'); hold on;
% 
% % Add identity line
% plot(theoretical, theoretical, 'r--', 'LineWidth', 1.5);
% 
% xlabel('Theoretical Log-normal Quantiles');
% ylabel('Empirical Quantiles');
% title('Q-Q Plot for Log-normal Fit');
% 
% grid on;


% disp("ANOVA test comparing models:");
% anova_results = anova(M1, M2, 'summary');
% disp(anova_results);

% % Extract log-likelihoods
% LL1 = M1.LogLikelihood;
% LL2 = M2.LogLikelihood;
% 
% % LR statistic
% LR = 2 * (LL2 - LL1);
% 
% % Degrees of freedom = # of extra parameters in Model 2
% df = M2.NumCoefficients - M1.NumCoefficients;
% 
% % p-value
% p = 1 - chi2cdf(LR, df);
% 
% fprintf('\nLikelihood Ratio Test:\n');
% fprintf('LR = %.3f, df = %d, p = %.5f\n', LR, df, p);


% %%
% figure;
% gscatter(T.orientation, T.certainty, T.uncertaintyLevel);
% xlabel('Orientation');
% ylabel('Confidence');
% title('Confidence across Orientation and Uncertainty');



%% Stats split by orientation

figure(9)

% orientations               = orientations;
n_uncertainty_levels       = uncertainty_levels;
n_orientations             = numel(orientations);
% resp_err_all               = resp_err_all;
% confidence_report_all      = confidence_report_all;

for i=1:n_uncertainty_levels
    for j=1:n_orientations

        errs_         = resp_err_all(i, j, :);
        conf_reports_ = confidence_report_all(i, j, :);
        
        countHC_ = sum(conf_reports_ == 1);
        countLC_ = sum(conf_reports_ == 0);
        
        subplot(n_uncertainty_levels, n_orientations, n_orientations*(i - 1) + j)
        bar([0, 1], [countHC_, countLC_])
        xticks([0 1])
        xticklabels(["HC", "LC"]);
        % ylim([0 15])

        if j == 1
            ylabel("count")
        end

        if i == 1
            title(orientations(j))
        end
    end
end

figure(10)

for i=1:n_uncertainty_levels
    for j=1:n_orientations

        errs_         = resp_err_all(i, j, :);
        conf_reports_ = confidence_report_all(i, j, :);
        
        stdHC_ = std( errs_(conf_reports_ == 1) );
        stdLC_ = std( errs_(conf_reports_ == 0) );

        % stdHC_ = mean( errs_(conf_reports_ == 1) );
        % stdLC_ = mean( errs_(conf_reports_ == 0) );
        
        subplot(n_uncertainty_levels, n_orientations, n_orientations*(i - 1) + j)
        bar([0, 1], [stdHC_, stdLC_])
        xticks([0 1])
        xticklabels(["HC", "LC"]);
        % ylim([0 60])
        
        if j == 1
            ylabel("std (data)")
        end
        
        if i == 1
            title(orientations(j))
        end
    end
end


%% Satatistical test of internal noise fluctuations HC and LC err split

errDiff_LCHC = nan + zeros(1, n_uncertainty_levels*n_orientations);
errLCs = nan + zeros(1, n_uncertainty_levels*n_orientations);
errHCs = nan + zeros(1, n_uncertainty_levels*n_orientations);

for i=1:n_uncertainty_levels
    for j=1:n_orientations

        errs_         = resp_err_all(i, j, :);
        conf_reports_ = confidence_report_all(i, j, :);
        
        errHC_ = errs_(conf_reports_ == 1);
        errLC_ = errs_(conf_reports_ == 0);
        
        if numel(errHC_) > 0 && numel(errLC_) > 0

            % For each subsample equal number of values and compute the
            % average error over many permutation

            if numel(errHC_) < numel(errLC_)
                y = mean( abs(errHC_) );

                x = 0;
                nitr = 100;
                for itr=1:nitr
                    x_samples = randsample(errLC_(:), numel(errHC_));
                    x = x + mean( abs(x_samples) ) ;
                end
                x = x / nitr;

            elseif numel(errHC_) > numel(errLC_)
                 x = mean( abs(errLC_) );

                 y = 0;
                 nitr = 100;
                 for itr=1:nitr
                    y_samples = randsample(errHC_(:), numel(errLC_));
                    y = y + mean( abs(y_samples) ) ;
                 end
                 y = y / nitr;
            else

                x = mean( abs(errLC_) );
                y = mean( abs(errHC_) );
            end

%             x = mean( abs(errLC_) );
%             y = mean( abs(errHC_) );

            errDiff = x - y;

            errDiff_LCHC( (i - 1)*n_uncertainty_levels + j ) = errDiff;
            errLCs( (i - 1)*n_uncertainty_levels + j )       = x;
            errHCs( (i - 1)*n_uncertainty_levels + j )       = y;
        else
            errDiff = nan;
        end

    end
end

[h,p,ci,stats] = ttest(errDiff_LCHC, 0)

figure
subplot(2, 2, 1)
hold on
histogram(errDiff_LCHC) %'BinEdges',-40:5:30
xline(0, 'LineStyle',"--")
xlabel("Abs error diff bw LC and HC")
ylabel("Count")
ylim([0 25])

subplot(2, 2, 2)
hold on
histogram(errLCs, DisplayName="LC") %'BinEdges',0:5:40, 
histogram(errHCs, DisplayName="HC") % 'BinEdges',0:5:40,
xline(0, 'LineStyle',"--")
xlabel("Abs error")
ylabel("Count")
legend
