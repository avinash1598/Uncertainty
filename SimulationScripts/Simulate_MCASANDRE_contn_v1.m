clear all
close all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')

orientations     = linspace(0, 180, 18); % linspace(0, 180, 18);
ntrials_per_ori  = 100;
sigma_e_levels   = linspace(0.01, 2, 10);  % Varying external noise or contrast
sigma_i          = 1;
varGain          = 0.5;
sigma_m          = 0.2;
Cc               = 0.5;

% Preallocate arrays
n_theta                  = length(orientations);
n_levels                 = length(sigma_e_levels);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(n_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(n_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(n_levels, n_theta, ntrials_per_ori);

% Simulation loop
for level = 1:numel(sigma_e_levels)
    sigma_e = sigma_e_levels(level);  % Set perceptual noise for this contrast

    for i = 1:numel(orientations)
        theta_true = orientations(i);   % True orientation
        trials = ntrials_per_ori;

        % Step1: Using estimate of this uncertainty, the subject estimates
        % the orientation
        % Internal estimate (Gaussian noise) - Note: this is not wraped
        % In actual behavioral data this will be wraped
        % theta_est = theta_true + sigma_p * randn(trials, 1);
        % Find doubly stochastic theta
        gain = gamrnd(1./varGain, varGain, [trials 1]);
        sigma_i_modulated = gain*sigma_i;
        sigma_theta = sqrt(sigma_e.^2 + sigma_i_modulated.^2);
        theta_est = theta_true + sigma_theta .* randn(trials, 1);

        assert(numel(sigma_theta) == trials);

        % Step1: Subject first gets an estimate of its uncertainty
        % Subject’s estimate of their uncertainty (meta-uncertainty)
        mu_log = log(sigma_theta.^2 ./ sqrt(sigma_m.^2 + sigma_theta.^2));
        sigma_log = sqrt(log(1 + (sigma_m.^2 ./ sigma_theta.^2)));
        sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
        
        % Confidence variable
        Vc = 1 ./ sigma_hat;
        
        % Store
        theta_true_all(level, i, :)         = theta_true;
        theta_resp_all(level, i, :)         = theta_est;
        confidence_report_all(level, i, :)  = Vc > Cc;
    end
end
   

% Plot the performance curves
% resp_err_all = abs(theta_resp_all - theta_true_all);
resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all_reshaped = reshape(resp_err_all, numel(sigma_e_levels), []);
confidence_report_all_reshaped = reshape(confidence_report_all, numel(sigma_e_levels), []);

data.err = resp_err_all_reshaped;
data.confReport = confidence_report_all_reshaped;
save('modelContOriData.mat', "data")

% raw_resp_err_all = theta_resp_all - theta_true_all;
% raw_resp_err_all_reshaped = reshape(raw_resp_err_all, numel(sigma_e_levels), []);

% resp_err_all_reshaped = squeeze(resp_err_all(:,1,:)); % Data for one angle only
% confidence_report_all_reshaped = squeeze(confidence_report_all(:,1,:));
        
figure

subplot(2, 3, 1)

x = mean(resp_err_all_reshaped, 2);
y = std(resp_err_all_reshaped, 0, 2);

errorbar(sigma_e_levels, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("External noise (\sigma_e)")
ylabel("Error")

subplot(2, 3, 2)

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

hold on

errorbar(sigma_e_levels, ...
    x_HC, y_HC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

errorbar(sigma_e_levels, ...
    x_LC, y_LC, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="Low confidence");

xlabel("External noise (\sigma_e)")
ylabel("Error")
legend
hold off

x1 = resp_HC(1, :); valid_idx = ~isnan(x1); x1 = x1(valid_idx);
x2 = resp_LC(1, :); valid_idx = ~isnan(x2); x2 = x2(valid_idx);
x3 = resp_HC(numel(sigma_e_levels), :); valid_idx = ~isnan(x3); x3 = x3(valid_idx);
x4 = resp_LC(numel(sigma_e_levels), :); valid_idx = ~isnan(x4); x4 = x4(valid_idx);

subplot(2, 3, 3)

% Behavioral variability
plot(sigma_e_levels, y, 'LineWidth', 2);
xlabel("External noise (\sigma_e)")
ylabel("Behavioral variability")

subplot(2, 3, 4)

% Behavioral variability
plot(sigma_e_levels, y_HC, 'LineWidth', 2, DisplayName="High confidence");
hold on
plot(sigma_e_levels, y_LC, 'LineWidth', 2, DisplayName="Low confidence");
xlabel("External noise (\sigma_e)")
ylabel("Behavioral variability")
legend
hold off



% Do bootstrap itterations and get difference in HC and LC error for each
% iterations
HC_idx = confidence_report_all_reshaped == 1;
LC_idx = confidence_report_all_reshaped == 0;

resp_HC = resp_err_all_reshaped;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_all_reshaped;
resp_LC(~LC_idx) = NaN;

nBoot = 10000;
someDummyArray = zeros(nBoot, numel(sigma_e_levels));
xvals = zeros(nBoot, numel(sigma_e_levels));

for i=1:nBoot
    for j=1:numel(sigma_e_levels)
        x_ = resp_HC(j, :);
        y_ = resp_LC(j, :);
        x_ = x_(~isnan(x_));
        y_ = y_(~isnan(y_));
        
        HC_star = x_(randi(numel(x_), numel(x_), 1));
        LC_star = y_(randi(numel(y_), numel(y_), 1));
        
        diff = mean(LC_star) - mean(HC_star);
        someDummyArray(i, j) = diff;
        xvals(i, j) = sigma_e_levels(j);
    end
end

subplot(2, 3, 5)

scatter(xvals(:), someDummyArray(:))
% plot(xvals(:), someDummyArray(:), LineWidth=2)
xlabel("\sigma_e levels")
ylabel("Difference in perceptual error b\w LC and HC")

%%
[corr_val, pval] = corr(xvals(:), someDummyArray(:), 'Type', 'Spearman');

doBootstrapTrendTest(sigma_e_levels, resp_HC, resp_LC)