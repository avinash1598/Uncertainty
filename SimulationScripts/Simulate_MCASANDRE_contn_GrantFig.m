%%%
% TODO:
% 1. Simulate effect of dispersion
% 2. Von misses distribution
%%%
clear all
close all


orientations     = linspace(0, 180, 18); % linspace(0, 180, 18);
ntrials_per_ori  = 25;
sigma_e_levels   = linspace(5, 30, 4);  % Varying external noise or contrast - simulate dispersion using this
sigma_i          = 8;
varGain          = 0.5;
sigma_m          = 2;
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
        % TODO: Check with von misses distribution
        
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
   
sz = size(theta_true_all);
stimOri = reshape(theta_true_all, sz(1), []);
reportedOri = reshape(theta_resp_all, sz(1), []);
confReports = reshape(confidence_report_all, sz(1), []);

wrappedReportedOri = mod(reportedOri, 180); % Wrap all the reported angles between 0 and 180
rawError = stimOri - wrappedReportedOri;
rawOriError = mod(rawError + 90, 180) - 90;


%% Figure 1 - Raw errors by group - Estimation error (Human subject)
% Define visual properties
colors = [0.85 0.33 0.10; 0 0.45 0.74; 0.85 0.33 0.10; 0 0.45 0.74]; % colorblind-friendly red/blue
alphas = [1, 1, 0.3, 0.3];

% Create figure
figure('Color','w');
pbaspect([3 5 1]) 
hold on

% Plot each group
for i = 1:numel(sigma_e_levels)
    grpOriErr = rawOriError(i, :);
    meanErr = mean(grpOriErr(~isnan(grpOriErr)));
    xPts = i + 0.4*(rand(1, numel(grpOriErr)) - 0.5);

    scatter(xPts, grpOriErr, 60, ...
        'MarkerEdgeColor', 'w', ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerFaceAlpha', alphas(i), ...
        'MarkerEdgeAlpha', alphas(i), ...
        'LineWidth', 0.5);

    plot([i - 0.35, i + 0.35], [meanErr, meanErr], LineStyle="-", LineWidth=2, Color='black')
end

% Labels and styling
ylabel("Estimation error (°)", 'FontSize', 14)
xlabel("Uncertainty level", 'FontSize', 14)
xticks(1:length(sigma_e_levels));
xticklabels(1:length(sigma_e_levels));
yticks([-90, 0, 90])
ylim([-90, 90])
xlim([0.3, 5])

% Improve figure aesthetics
set(gca, 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out', 'Box', 'off')
hold off

% Better to use this with alphadata
exportgraphics(gcf, 'Model1.eps', 'ContentType', 'vector')
% print(gcf, 'Model1', '-depsc', '-r300')  % -r300 = 300 DPI


%% Figure 2 - Error by uncertainty - take MSE error instead (maybe)
% === Figure 2 ===
% Clean data
rawErr = rawOriError(:);
errHC = rawErr(confReports == 1);
errLC = rawErr(confReports == 0);
errHC = errHC(~isnan(errHC));
errLC = errLC(~isnan(errLC));

abErrHC = abs(errHC);
abErrLC = abs(errLC);


% Prepare values
means = [mean(abErrLC), mean(abErrHC)];
sems  = [std(abErrLC)/sqrt(numel(abErrLC)), std(abErrHC)/sqrt(numel(abErrHC))]; % No bias in this plot needs to be shown
colors = [0.85 0.33 0.10; 0.47 0.67 0.19]; 

% Create figure
figure('Color','w', 'Units','inches', 'Position',[1, 1, 4.5, 5]);
pbaspect([3 5 1]) 
hold on

labels = ["Uncertain", "Certain"];
% Bar plot with custom colors
for i = 1:2
    b = bar(i+0.5, means(i), 0.6, DisplayName=labels(i)); % 0.6 = bar width
    b.FaceColor = colors(i,:);
    b.EdgeColor = 'none';
end

% Add error bars
errorbar([1.5 2.5], means, sems, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10, HandleVisibility='off');

% Styling
ylabel("|Estimation Error| (°)", 'FontSize', 14)
yticks([0, 30])
xticks([])
% xticklabels({'Low Conf.', 'High Conf.'})
lgd = legend('\color[rgb]{0.85, 0.33, 0.10} Uncertain','\color[rgb]{0.47, 0.67, 0.19} Certain');
lgd.Box = "off";
set(gca, 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out', 'Box', 'off')
% ylim([0, max(means + sems) * 1.2])  % Leave headroom
ylim([0, 30])  % Leave headroom

% Grid and finish
hold off

exportgraphics(gcf, 'Model2.eps', 'ContentType', 'vector')
% print(gcf, 'Model2', '-depsc', '-r300')  % -r300 = 300 DPI


