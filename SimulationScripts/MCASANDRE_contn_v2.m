clear all
close all

sigma_reward      = 0.5;
reward_K          = 0.2;  % Constant reward for low confidence
sigma_e           = linspace(0.01, 2, 100);

modelParams.sigma_e                   = 0.7;
modelParams.sigma_i                   = 1;
modelParams.varGain                   = 0.5;
modelParams.internalNoiseSamplesCnt   = 200;
modelParams.Cc                        = 0.5;
modelParams.sigma_m                   = 0.2;

jointPDF_HC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
jointPDF_LC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
sigma_theta = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);

for idx=1:numel(sigma_e)
    modelParams.sigma_e = sigma_e(idx);     
    retData = getLLhChoice_MCASANDRE_Contn(modelParams);

    jointPDF_HC(idx, :) =  retData.pHC;
    jointPDF_LC(idx, :) =  retData.pLC;
    sigma_theta(idx, :) =  retData.sigma_theta; % sigma_theta is same for all
end

% Reward computation
max_sigma = max(max(sigma_theta));
x = linspace( - 4*max_sigma, 4*max_sigma, 1000);

N = size(sigma_theta, 1);
M = size(sigma_theta, 2);
K = length(x);
x_exp = reshape(x, 1, 1, K);  % (1×1×K)
sigma_exp = reshape(sigma_theta, N, M, 1);  % (N×M×1)

rewardFn     = normpdf(x, 0, sigma_reward);
rewardFn_exp = reshape(rewardFn, 1, 1, []);
errorFn      = normpdf(x_exp, 0, sigma_exp);

expRewardHC  = sum(rewardFn_exp.*errorFn, 3);
expRewardLC  = sum(reward_K.*errorFn, 3);
netExpReward = jointPDF_HC.*expRewardHC + jointPDF_LC.*expRewardLC;

figure

subplot(2, 2, 4)
hold on

x1 = sigma_theta(:);
y1 = expRewardHC(:);
y2 = expRewardLC(:);
y3 = netExpReward(:);

[x1_sorted, idx] = sort(x1);
y1_sorted        = y1(idx);
y2_sorted        = y2(idx);
y3_sorted        = y3(idx);

plot(x1_sorted, y1_sorted, DisplayName='HC reward', LineWidth=2)
plot(x1_sorted, y2_sorted, DisplayName='LC reward', LineWidth=2)
plot(x1_sorted, y3_sorted, DisplayName='Net reward', LineWidth=3, LineStyle='--')
% scatter(x1_sorted, y1_sorted, DisplayName='HC reward')
% scatter(x1_sorted, y2_sorted, DisplayName='LC reward')

xlabel("\sigma_theta")
ylabel("Expected reward")
legend

hold off


subplot(2, 2, 1)

shapeParam = 1./modelParams.varGain;
scaleParam = modelParams.varGain;
sigma_i = modelParams.sigma_i*gamrnd(shapeParam, scaleParam, [1 100]); % Gain vector for all trials

x = 1:numel(sigma_e);
hold on
plot(x, sigma_e, DisplayName='External noise', LineWidth=1.5);

x1 = repmat(x, [numel(sigma_i) 1]);
y1 = repmat(sigma_i, [numel(sigma_e) 1])';
mean_y1 = mean(y1, 1);

x1 = reshape(x1, [], 1);
y1 = reshape(y1, [], 1);
scatter(x1, y1, HandleVisibility='off', ...
    DisplayName='Internal noise', ...
    MarkerFaceAlpha = 0.1, ...
    MarkerEdgeAlpha = 0.1);

scatter(x, mean_y1, 10, "filled", "o", ...
    MarkerEdgeColor='b', MarkerFaceColor='b', ...
    HandleVisibility='on', DisplayName='Internal noise');

ylabel("Noise")

legend

hold off


subplot(2, 2, 2)

mean_sigma_theta  = mean(sigma_theta, 2);
std_sigma_theta   = std(sigma_theta, 0, 2);

hold on

fill([sigma_e'; flipud(sigma_e')], ...
     [mean_sigma_theta + std_sigma_theta; flipud(mean_sigma_theta - std_sigma_theta)], ...
     [1 0.8 0.8], ...           % light red color
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4, ...
     'HandleVisibility','off', ...
     'DisplayName', 'LC ± std');

% Plot means
plot(sigma_e', mean_sigma_theta, 'b-', 'LineWidth', 1.5);

xlabel("External noise")
ylabel("Expected [Sigma theta]")

hold off


subplot(2, 2, 3)

% Maybe split sigma_theta by HC and LC 
% calculate expected error: pHC*sigma_theta

% Calculate expected standard deviation for low confidence and high
% confidence choices
exp_sig_theta_HC = sum(sigma_theta.*jointPDF_HC, 2) ./ sum(jointPDF_HC, 2);
exp_sig_theta_LC = sum(sigma_theta.*jointPDF_LC, 2) ./ sum(jointPDF_LC, 2);

% Varaince
exp_sig_theta2_HC = sum((sigma_theta.^2) .* jointPDF_HC, 2) ./ sum(jointPDF_HC, 2);
std_sig_theta_HC = sqrt( exp_sig_theta2_HC - exp_sig_theta_HC.^2 );

exp_sig_theta2_LC = sum((sigma_theta.^2) .* jointPDF_LC, 2) ./ sum(jointPDF_LC, 2);
std_sig_theta_LC = sqrt( exp_sig_theta2_LC - exp_sig_theta_LC.^2 );


hold on

% LC shaded area
fill([sigma_e'; flipud(sigma_e')], ...
     [exp_sig_theta_LC + std_sig_theta_LC; flipud(exp_sig_theta_LC - std_sig_theta_LC)], ...
     [1 0.8 0.8], ...           % light red color
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4, ...
     'HandleVisibility','off', ...
     'DisplayName', 'LC ± std');

% HC shaded area
fill([sigma_e'; flipud(sigma_e')], ...
     [exp_sig_theta_HC + std_sig_theta_HC; flipud(exp_sig_theta_HC - std_sig_theta_HC)], ...
     [0.8 0.8 1], ...           % light blue color
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4, ...
     'HandleVisibility','off',...
     'DisplayName', 'HC ± std');

% Plot means
plot(sigma_e', exp_sig_theta_HC, 'b-', 'LineWidth', 1.5, 'DisplayName', 'HC');
plot(sigma_e', exp_sig_theta_LC, 'r-', 'LineWidth', 1.5, 'DisplayName', 'LC');

xlabel("External noise")
ylabel("Expected [Sigma theta]")
legend

hold off






