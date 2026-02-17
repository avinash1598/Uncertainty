%%%
% TODO:
% 1. Simulate effect of dispersion
% 2. Von misses distribution

% Run following experiments:
% Need better estimate of sigma_i and varGain
% Does Cc depend upon Reward

%%%
clear all
close all


orientation      = 0;
ntrials_per_ori  = 500;
sigma_e          = 9;  
sigma_i          = 14;
varGain          = 0.5;
sigma_m          = 2;
% Cc               = 0.5;
Ccs              = linspace(0, 0.5, 100);

% Only record data which would actually be recorded during experiment
theta_true_all        = zeros(1, ntrials_per_ori);
theta_resp_all        = zeros(1, ntrials_per_ori); % Recorded theta based on user response
confidence_report_all = zeros(1, ntrials_per_ori);

% Simulation loop
theta_true = orientation;   % True orientation
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
theta_true_all(:)         = theta_true;
theta_resp_all(:)         = theta_est;
confidence_report_all(:)  = Vc > Ccs(10);

HC_tot_reward_List         = zeros(size(Ccs));
LC_tot_reward_List         = zeros(size(Ccs));
total_reward_all_LC        = zeros(size(Ccs));
total_reward_all_HC        = zeros(size(Ccs));
total_reward_random_guess  = zeros(size(Ccs));

for i=1:numel(Ccs)
    Cc = Ccs(i);
    confidence_report_all(:)  = Vc > Cc;
    rewardVals = calcReward(theta_true_all, theta_resp_all, confidence_report_all);

    totalReward_HC = sum(rewardVals(confidence_report_all == 1));
    totalReward_LC = sum(rewardVals(confidence_report_all == 0));

    HC_tot_reward_List(i) = totalReward_HC;
    LC_tot_reward_List(i) = totalReward_LC;

    % All HC
    conf_report = 1 + zeros(size(confidence_report_all));
    rewardVals = calcReward(theta_true_all, theta_resp_all, conf_report);
    total_reward_all_HC(i) = sum(rewardVals);

    % All LC
    conf_report = zeros(size(confidence_report_all));
    rewardVals = calcReward(theta_true_all, theta_resp_all, conf_report);
    total_reward_all_LC(i) = sum(rewardVals);
    
    % Random guesses
    randVals = rand(size(confidence_report_all));
    conf_report = randVals > 0.5;
    rewardVals = calcReward(theta_true_all, theta_resp_all, conf_report);
    total_reward_random_guess(i) = sum(rewardVals);

end

figure
hold on
plot(Ccs, LC_tot_reward_List, DisplayName="LC", LineWidth=2)
plot(Ccs, HC_tot_reward_List, DisplayName="HC", LineWidth=2)
plot(Ccs, total_reward_all_HC, DisplayName="all HC", LineWidth=2)
plot(Ccs, total_reward_all_LC, DisplayName="all LC", LineWidth=2)
plot(Ccs, total_reward_random_guess, DisplayName="Random guess", LineWidth=2)
xlabel("Cc")
ylabel("Total reward")
legend
hold off

% rewardVals = calcReward(theta_true_all, theta_resp_all, confidence_report_all);
% 
% err_all = theta_true_all - theta_resp_all;
% 
% figure
% hold on
% scatter(err_all(confidence_report_all == 1), rewardVals(confidence_report_all == 1), DisplayName="HC")
% scatter(err_all(confidence_report_all == 0), rewardVals(confidence_report_all == 0), DisplayName="LC")
% xlabel("Error")
% ylabel("Reward")
% legend
% yline(0, LineStyle="--", HandleVisibility="off")
% % xlim([-20, 20])
% hold off

% for i=1:numel(Ccs)
%     Cc = Ccs(i);
%     confReport = Vc > Cc;
% end



% sz = size(theta_true_all);
% stimOri = reshape(theta_true_all, sz(1), []);
% reportedOri = reshape(theta_resp_all, sz(1), []);
% confReports = reshape(confidence_report_all, sz(1), []);
% 
% wrappedReportedOri = mod(reportedOri, 180); % Wrap all the reported angles between 0 and 180
% rawError = stimOri - wrappedReportedOri;
% rawOriError = mod(rawError + 90, 180) - 90;

% p = linspace(0, 1, 20);
% mu = 
% sigma = sigma_m;
% x = logninv(p,mu,sigma)

% Calculate reward
function reward = calcReward(trueOri, reportedOri, confReport)

% Note: reported orientation is already pi-periodic, true orientation as
% well
maxTolerableError = 25; % In degrees
absPerceptualError = abs(trueOri - reportedOri);
absPerceptualError = min(absPerceptualError, 180 - absPerceptualError);
reward = zeros(size(absPerceptualError));

% y1 = 9;
% x1 = 7.28;
% c1 = 1;
% c2 = 0.3;
% m1 = 1/y1;
% m2 = m1 - (c1 - c2)/x1;

y1 = 12;
x1 = 10;
c1 = 1;
c2 = 0.3;
m1 = c1/y1;
m2 = c2 / maxTolerableError; %m1 - (c1 - c2)/x1;

% High confidence
gHC = - m1 * absPerceptualError + c1;
reward(:) = gHC; %6*gHC
idx = (absPerceptualError <= y1) & (confReport == 1);
reward(idx) = gHC(idx);
idx = (absPerceptualError > y1) & (confReport == 1);
reward(idx) = 0*gHC(idx); %0.3*gHC(idx);

% Low confidence
gLC = - m2 * absPerceptualError + c2;
reward(confReport == 0) = gLC(confReport == 0);

% HC error
idx = ( absPerceptualError > maxTolerableError ) & (confReport == 1);
reward(idx) = -0.8; %-0.4; 

% LC error
idx = ( absPerceptualError > maxTolerableError ) & (confReport == 0);
reward(idx) = 0;

end