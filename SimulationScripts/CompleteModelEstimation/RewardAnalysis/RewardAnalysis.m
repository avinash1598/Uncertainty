% TODO: have a utility function that pulls this out

close all
clear all

addpath('LL_scripts/')

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/data/COR15.mat'); % Changed reward function - affects perceptual variance - more relaxed
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/data/COR31.mat'); % Changed reward function - affects perceptual variance - more relaxed
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR33.mat');

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


%%
confReports = confidence_report_all(:);
respErr = resp_err_all(:);
rvOriErr     = -90:3:90;

respErrHC = respErr(confReports == 1);
respErrLC = respErr(confReports == 0);

[pdf, edges] = histcounts(respErr, ...
        'Normalization', 'pdf', ...
        'BinEdges', rvOriErr);

% Compute CDF
centers = edges(1:end-1) + diff(edges)/2;
cdf = cumsum(pdf .* diff(edges));
cdf = cdf / cdf(end); % Normalize to ensure it ends at 1 (numerical precision)

[cdf_75, idx_75] = min(abs(cdf - 0.75));
[cdf_25, idx_25] = min(abs(cdf - 0.25));

figure
subplot(2, 3, 1)
hold on
histogram(respErrHC, BinEdges=rvOriErr, DisplayName="HC", Normalization="pdf", FaceAlpha=0.5)
histogram(respErrLC, BinEdges=rvOriErr, DisplayName="LC", Normalization="pdf", FaceAlpha=0.5)
legend
hold off

subplot(2, 3, 2)
histogram(respErrHC, BinEdges=rvOriErr, DisplayName="HC", Normalization="pdf", FaceAlpha=0.5)
title("HC")

subplot(2, 3, 3)
histogram(respErrLC, BinEdges=rvOriErr, DisplayName="LC", Normalization="pdf", FaceAlpha=0.5)
title("LC")

subplot(2, 3, 4)
hold on
histogram(respErr, BinEdges=rvOriErr, Normalization="pdf", FaceAlpha=0.5)
xline(centers(idx_25), LineStyle="--", LineWidth=1.5)
xline(centers(idx_75), LineStyle="--", LineWidth=1.5)
title("All")
hold off

disp(centers(idx_25))
disp(centers(idx_75))

% 50% performance tolerance
centers(idx_75) - centers(idx_25)

% ideal range -9, 9 => this should be the intersection point b/w LC and HC
% reward function

subplot(2, 3, 5)
bar([1, 2], [sum(confidence_report_all(:) == 1) sum(confidence_report_all(:) == 0) ])
xticklabels(["HC", "LC"])
ylabel("trials count")

%%
%%%
% TODO:
% 1. Simulate effect of dispersion
% 2. Von misses distribution

% Run following experiments:
% Need better estimate of sigma_i and varGain
% Does Cc depend upon Reward

%%%
clear all
% close all

orientations     = 0:10:180; % linspace(0, 180, 18);
ntrials_per_ori  = 25;
sigma_s          = [12.1392   13.9612   19.8033   19.6392   22.2473   26.4453   23.1557   26.5658];
scale            = 358;
sigma_meta       = 5.2349;
Cc               = 0.06; 

% Preallocate arrays
n_theta                  = numel(orientations);
uncertainty_levels       = numel(sigma_s);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all        = zeros(uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
VcEstimates           = zeros(uncertainty_levels, n_theta, ntrials_per_ori);
theta_resp_all_RG     = zeros(uncertainty_levels, n_theta, ntrials_per_ori);

for l=1:uncertainty_levels
    for i = 1:n_theta
        theta_true = orientations(i);   % True orientation
        trials = ntrials_per_ori;
        
        % Multiplicative
        shape = sigma_s(l).^2 / scale; % divide by scale so that mean is sigma_s
        gain = gamrnd(shape, scale, [trials 1]);
        sigma_m = sqrt( gain );
        mean_m = theta_true;

        % Wrap the angle at the plotting stage. Note: warapping should be
        % performed according to the true angle.
        theta_est = mean_m + sigma_m .* randn(trials, 1);
        % Since this is orientation, wrap the angle between 0 and 180
        theta_est = mod(theta_est, 180); 
        % Are the actual distribution near 0 and 180 captured by this in simulation
    
        assert(numel(sigma_m) == trials);
        
        % Step1: Subject first gets an estimate of its uncertainty
        % Subject’s estimate of their uncertainty (meta-uncertainty)
        mu_log = log(sigma_m.^2 ./ sqrt(sigma_meta.^2 + sigma_m.^2));
        sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m.^2)));
        sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
        
        % Confidence variable
        Vc = 1 ./ sigma_hat;
        
        % Store
        theta_true_all(l, i, :)         = theta_true;
        theta_resp_all(l, i, :)         = theta_est;
        VcEstimates(l, i, :)            = Vc;
        confidence_report_all(l, i, :)  = Vc > Cc;
        theta_resp_all_RG(l, i, :)      = unifrnd(0, 179, [trials, 1]);
    end
end

%% Plot results
% Plot the performance curves
resp_err_all = (theta_resp_all - theta_true_all);
resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error

resp_err_all_reshaped = reshape(resp_err_all, uncertainty_levels, []);
confidence_report_all_reshaped = reshape(confidence_report_all, uncertainty_levels, []);

figure

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
errorbar(sigma_s, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("\sigma_s(s)")
ylabel("Error")

subplot(2, 3, 2)

% Behavioral variability
scatter(sigma_s, y, "filled");
hold on
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim, LineWidth=1.5);
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
hold off

subplot(2, 3, 3)

% Behavioral variability
scatter(sigma_s, y_HC, "filled", DisplayName="High confidence");
hold on
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_HC, LineWidth=1.5, HandleVisibility="off");
scatter(sigma_s, y_LC, "filled", DisplayName="Low confidence");
% plot(mean_sigma_s_stim, anlytcl_sigma_m_stim_LC, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (measurement noise)")
ylabel("\sigma_m(s) (sensory noise)")
legend
hold off


%% Simulation
confReports = confidence_report_all(:);
respErr = resp_err_all(:);
rvOriErr     = -90:3:90;

respErrHC = respErr(confReports == 1);
respErrLC = respErr(confReports == 0);

[pdf, edges] = histcounts(respErr, ...
        'Normalization', 'pdf', ...
        'BinEdges', rvOriErr);

% Compute CDF
centers = edges(1:end-1) + diff(edges)/2;
cdf = cumsum(pdf .* diff(edges));
cdf = cdf / cdf(end); % Normalize to ensure it ends at 1 (numerical precision)

[cdf_75, idx_75] = min(abs(cdf - 0.75));
[cdf_25, idx_25] = min(abs(cdf - 0.25));

figure
subplot(2, 3, 1)
hold on
histogram(respErrHC, BinEdges=rvOriErr, DisplayName="HC", Normalization="pdf", FaceAlpha=0.5)
histogram(respErrLC, BinEdges=rvOriErr, DisplayName="LC", Normalization="pdf", FaceAlpha=0.5)
legend
hold off

subplot(2, 3, 2)
histogram(respErrHC, BinEdges=rvOriErr, DisplayName="HC", Normalization="pdf", FaceAlpha=0.5)
title("HC")

subplot(2, 3, 3)
histogram(respErrLC, BinEdges=rvOriErr, DisplayName="LC", Normalization="pdf", FaceAlpha=0.5)
title("LC")

subplot(2, 3, 4)
hold on
histogram(respErr, BinEdges=rvOriErr, Normalization="pdf", FaceAlpha=0.5)
xline(centers(idx_25), LineStyle="--", LineWidth=1.5)
xline(centers(idx_75), LineStyle="--", LineWidth=1.5)
title("All")
hold off

%%
resp_err_all_reshaped            = resp_err_all(:);
confidence_report_all_reshaped   = confidence_report_all(:);
VcEstimates_reshaped             = VcEstimates(:);

Ccs                        = linspace(0, 0.5, 1000);
HC_tot_reward_List         = zeros(size(Ccs));
LC_tot_reward_List         = zeros(size(Ccs));
tot_reward_HC_LC           = zeros(size(Ccs));
total_reward_all_LC        = zeros(size(Ccs));
total_reward_all_HC        = zeros(size(Ccs));
total_reward_all_LC_RG     = zeros(size(Ccs));
total_reward_all_HC_RG      = zeros(size(Ccs));
total_reward_random_guess  = zeros(size(Ccs));

for i=1:numel(Ccs)
    Cc_ = Ccs(i);
    confidence_report_all_  = VcEstimates_reshaped > Cc_;
    rewardVals = calcReward(theta_true_all(:), theta_resp_all(:), confidence_report_all_(:));

    totalReward_HC = sum(rewardVals(confidence_report_all_ == 1));
    totalReward_LC = sum(rewardVals(confidence_report_all_ == 0));

    tot_reward_HC_LC(i) = totalReward_HC + totalReward_LC;
    HC_tot_reward_List(i) = totalReward_HC;
    LC_tot_reward_List(i) = totalReward_LC;

    % All HC
    conf_report = 1 + zeros(size(confidence_report_all_));
    rewardVals = calcReward(theta_true_all(:), theta_resp_all(:), conf_report);
    total_reward_all_HC(i) = sum(rewardVals);

    % All LC
    conf_report = zeros(size(confidence_report_all_));
    rewardVals = calcReward(theta_true_all(:), theta_resp_all(:), conf_report);
    total_reward_all_LC(i) = sum(rewardVals);
    
    % Random guesses
    randVals = rand(size(confidence_report_all_));
    conf_report = randVals > 0.5;
    rewardVals = calcReward(theta_true_all(:), theta_resp_all(:), conf_report);
    total_reward_random_guess(i) = sum(rewardVals);

    % All LC but random guesses
    conf_report = 1 + zeros(size(confidence_report_all_));
    rewardVals = calcReward(theta_true_all(:), theta_resp_all_RG(:), conf_report);
    total_reward_all_LC_RG(i) = sum(rewardVals);

    % All HC but random guesses
    conf_report = zeros(size(confidence_report_all_));
    rewardVals = calcReward(theta_true_all(:), theta_resp_all_RG(:), conf_report);
    total_reward_all_HC_RG(i) = sum(rewardVals);

end

figure
hold on
plot(Ccs, tot_reward_HC_LC, DisplayName="Optimal", LineWidth=2)
plot(Ccs, LC_tot_reward_List, DisplayName="LC optimal", LineWidth=2, lineStyle="--")
plot(Ccs, HC_tot_reward_List, DisplayName="HC optimal", LineWidth=2, lineStyle="--")
plot(Ccs, total_reward_all_HC, DisplayName="all HC", LineWidth=2)
plot(Ccs, total_reward_all_LC, DisplayName="all LC", LineWidth=2)
plot(Ccs, total_reward_random_guess, DisplayName="Random guess", LineWidth=2)
% plot(Ccs, total_reward_all_LC_RG, DisplayName="All LC (Random guess)", LineWidth=2)
% plot(Ccs, total_reward_all_HC_RG, DisplayName="All HC (Random guess)", LineWidth=2)
xlabel("Cc")
ylabel("Total reward")
legend
hold off


% Step1: choose LC function i.e. total reward possible for low confidence
% Step2: tune HC function

% Will these parameters lead to 50:50 split bw HC and LC?

% Calculate reward
function reward = calcReward(trueOri, reportedOri, confReport)

% Note: reported orientation is already pi-periodic, true orientation as
% well
maxTolerableError = 30; %30; % In degrees
% absPerceptualError = abs(trueOri - reportedOri);
% absPerceptualError = min(absPerceptualError, 180 - absPerceptualError);
rawError = reportedOri - trueOri;
rawOriError = mod(rawError + 90, 180) - 90;
absPerceptualError = abs(rawOriError);

reward = zeros(size(absPerceptualError));

% y1 = 12 + (14 - 12)*abs(sind(2*trueOri)); % 12;
% x1 = 11 + (13 - 11)*abs(sind(2*trueOri)); % ideally 9 but its okay
y1 = 18; %18; %15; % 12;
x1 = 16; %16.1; %14; % ideally 9 but its okay
c1 = 0.15; %0.5; %0.2
m1 = - c1./y1;
m2 =  ( m1.*x1 + c1 ) ./ (x1 - maxTolerableError); % c2 / maxTolerableError; %m1 - (c1 - c2)/x1;
c2 = - m2.*maxTolerableError;

% High confidence
gHC = m1 .* absPerceptualError + c1;
reward(:) = gHC;
idx = (absPerceptualError <= y1) & (confReport == 1);
reward(idx) = gHC(idx);
idx = (absPerceptualError > y1) & (confReport == 1);
reward(idx) = 1.2*gHC(idx); %0.3*gHC(idx);

% Low confidence
gLC = m2 .* absPerceptualError + c2;
reward(confReport == 0) = gLC(confReport == 0);

% HC error
idx = ( absPerceptualError > maxTolerableError ) & (confReport == 1);
reward(idx) = -1.5*c1; %-0.3*c1; %-1.6*c1; 

% LC error
idx = ( absPerceptualError > maxTolerableError ) & (confReport == 0);
reward(idx) = 0;

% reward = 0.1*reward;

end




% Tese parameters look nice!
% % Calculate reward
% function reward = calcReward(trueOri, reportedOri, confReport)
% 
% % Note: reported orientation is already pi-periodic, true orientation as
% % well
% maxTolerableError = 25; % In degrees
% % maxTolerableErrorHC = 20;
% maxTolerableErrorHC = maxTolerableError;
% absPerceptualError = abs(trueOri - reportedOri);
% absPerceptualError = min(absPerceptualError, 180 - absPerceptualError);
% reward = zeros(size(absPerceptualError));
% 
% y1 = 12;
% x1 = 11; % ideally 9 but its okay
% c1 = 0.2;
% m1 = - c1/y1;
% m2 =  ( m1*x1 + c1 ) / (x1 - maxTolerableError); % c2 / maxTolerableError; %m1 - (c1 - c2)/x1;
% c2 = - m2*maxTolerableError;
% 
% % High confidence
% gHC = m1 * absPerceptualError + c1;
% reward(:) = gHC;
% idx = (absPerceptualError <= y1) & (confReport == 1);
% reward(idx) = gHC(idx);
% idx = (absPerceptualError > y1) & (confReport == 1);
% reward(idx) = 0.5*gHC(idx);
% 
% % Low confidence
% gLC = m2 * absPerceptualError + c2;
% reward(confReport == 0) = gLC(confReport == 0);
% 
% % HC error
% % idx = ( absPerceptualError > maxTolerableError ) & (confReport == 1);
% % reward(idx) = -0.2; %-0.4; 
% idx = ( absPerceptualError > maxTolerableErrorHC ) & (confReport == 1);
% reward(idx) = -1.6*c1; %-0.2; %-0.4; 
% 
% % LC error
% idx = ( absPerceptualError > maxTolerableError ) & (confReport == 0);
% reward(idx) = 0;
% 
% % reward = 0.66*reward;
% end