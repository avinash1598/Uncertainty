clear all
close all

sigma_reward      = 0.5;
reward_K          = 0.2;  % Constant reward for low confidence
sigma_e           = linspace(0.01, 2, 100);

modelParams.sigma_e                   = 0.7;
modelParams.sigma_i                   = 1;   % 1
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

subplot(2, 2, 1)

shapeParam = 1./modelParams.varGain;
scaleParam = modelParams.varGain;
sigma_i_ = 2*sigma_e.*gamrnd(shapeParam, scaleParam, [101 1]); % Gain vector for all trials

x = 1:numel(sigma_e);
hold on
plot(x, sigma_e, DisplayName='External noise', LineWidth=1.5);

x1_ = repmat(x, [size(sigma_i_, 1) 1]);
y1_ = sigma_i_;
mean_y1_ = mean(y1_, 1);
% 
% x1 = reshape(x1, [], 1);
% y1 = reshape(y1, [], 1);
scatter(x1_(:), y1_(:), HandleVisibility='off', ...
    DisplayName='Internal noise', ...
    MarkerFaceAlpha = 0.1, ...
    MarkerEdgeAlpha = 0.1);

scatter(x, mean_y1_, 10, "filled", "o", ...
    MarkerEdgeColor='b', MarkerFaceColor='b', ...
    HandleVisibility='on', DisplayName='Internal noise');

ylabel("Noise")

legend

hold off

subplot(2, 2, 2)
% Negative covariation
sigma_i0 = 0.9;
k = 0.001;

shapeParam = 1./modelParams.varGain;
scaleParam = modelParams.varGain;
sigma_i_ = sqrt( k ./ (sigma_e.^2) + sigma_i0^2 );
sigma_i_ = sigma_i_.*gamrnd(shapeParam, scaleParam, [101 1]); % Gain vector for all trials

x = 1:numel(sigma_e);
hold on
plot(x, sigma_e, DisplayName='External noise', LineWidth=1.5);

x1_ = repmat(x, [size(sigma_i_, 1) 1]);
y1_ = sigma_i_;
mean_y1_ = mean(y1_, 1);

scatter(x1_(:), y1_(:), HandleVisibility='off', ...
    DisplayName='Internal noise', ...
    MarkerFaceAlpha = 0.1, ...
    MarkerEdgeAlpha = 0.1);

scatter(x, mean_y1_, 10, "filled", "o", ...
    MarkerEdgeColor='b', MarkerFaceColor='b', ...
    HandleVisibility='on', DisplayName='Internal noise');

ylabel("Noise")

legend

hold off







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
% sigma_i = 2*sigma_e.*gamrnd(shapeParam, scaleParam, [100 1]); % Gain vector for all trials

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







sigma_e           = linspace(0.01, 3, 100);
sigma_m           = linspace(0.2, 20, 2);

modelParams.sigma_e                   = 0.7;
modelParams.sigma_i                   = 1;
modelParams.varGain                   = 0.5;
modelParams.internalNoiseSamplesCnt   = 200;
modelParams.Cc                        = 1;
modelParams.sigma_m                   = 0.2;

figure

for midx=1:numel(sigma_m)
    jointPDF_HC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
    jointPDF_LC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
    sigma_theta = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);

    for idx=1:numel(sigma_e)
        modelParams.sigma_m = sigma_m(midx);
        modelParams.sigma_e = sigma_e(idx);     
        retData = getLLhChoice_MCASANDRE_Contn(modelParams);
    
        jointPDF_HC(idx, :) =  retData.pHC;
        jointPDF_LC(idx, :) =  retData.pLC;
        sigma_theta(idx, :) =  retData.sigma_theta; % sigma_theta is same for all
    end


    subplot(2, 1, 1)

    mean_sigma_theta  = mean(sigma_theta, 2);
    
    hold on

    plot(sigma_e', mean_sigma_theta, 'b-', 'LineWidth', 1.5, ...
        DisplayName=['\sigma_m = ' num2str(sigma_m(midx))]);
    
    xlabel("External noise")
    ylabel("Expected [Sigma theta]")
    legend

    hold off
    


    subplot(2, 1, 2)

    hold on

    exp_sig_theta_HC = sum(sigma_theta.*jointPDF_HC, 2) ./ sum(jointPDF_HC, 2);
    exp_sig_theta_LC = sum(sigma_theta.*jointPDF_LC, 2) ./ sum(jointPDF_LC, 2);
    
    % Plot means
    if midx > 1
        lineStyle="--";
    else
        lineStyle="-";
    end

    plot(sigma_e', exp_sig_theta_HC, 'b-', 'LineWidth', 1.5, 'DisplayName', 'HC', ...
        lineStyle=lineStyle, ...
        DisplayName=['HC \sigma_m = ' num2str(sigma_m(midx))]);
    plot(sigma_e', exp_sig_theta_LC, 'r-', 'LineWidth', 1.5, 'DisplayName', 'LC', ...
        lineStyle=lineStyle, ...
        DisplayName=['LC \sigma_m = ' num2str(sigma_m(midx))]);
    
    xlabel("External noise")
    ylabel("Expected [Sigma theta]")
    legend
    
    hold off

end







sigma_e           = linspace(0.01, 3, 100);

modelParams.sigma_e                   = 0.7;
modelParams.sigma_i                   = 1;
modelParams.varGain                   = 0.5;
modelParams.internalNoiseSamplesCnt   = 200;
modelParams.Cc                        = 1;
modelParams.sigma_m                   = 0.2;

figure
sgtitle("Internal and External noise covaries (internal << external)")

jointPDF_HC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
jointPDF_LC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
sigma_theta = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);

for idx=1:numel(sigma_e)
    modelParams.sigma_i = sqrt((0.2*sigma_e(idx)^2)); % sigma_i less than sigma_e always
    modelParams.sigma_e = sigma_e(idx);     
    retData = getLLhChoice_MCASANDRE_Contn(modelParams);
    
    jointPDF_HC(idx, :) =  retData.pHC;
    jointPDF_LC(idx, :) =  retData.pLC;
    sigma_theta(idx, :) =  retData.sigma_theta; % sigma_theta is same for all
end


subplot(2, 1, 1)

mean_sigma_theta  = mean(sigma_theta, 2);

hold on

plot(sigma_e', mean_sigma_theta, 'b-', 'LineWidth', 1.5);

xlabel("External noise")
ylabel("Expected [Sigma theta]")

hold off



subplot(2, 1, 2)

hold on

exp_sig_theta_HC = sum(sigma_theta.*jointPDF_HC, 2) ./ sum(jointPDF_HC, 2);
exp_sig_theta_LC = sum(sigma_theta.*jointPDF_LC, 2) ./ sum(jointPDF_LC, 2);

plot(sigma_e', exp_sig_theta_HC, 'b-', 'LineWidth', 1.5, 'DisplayName', 'HC', ...
    DisplayName='HC');
plot(sigma_e', exp_sig_theta_LC, 'r-', 'LineWidth', 1.5, 'DisplayName', 'LC', ...
    DisplayName='LC');

xlabel("External noise")
ylabel("Expected [Sigma theta]")
legend

hold off



sigma_e           = linspace(0.01, 3, 100);

modelParams.sigma_e                   = 0.7;
modelParams.sigma_i                   = 1;
modelParams.varGain                   = 0.5;
modelParams.internalNoiseSamplesCnt   = 200;
modelParams.Cc                        = 1;
modelParams.sigma_m                   = 0.2;

figure
sgtitle("Internal and External noise covaries (internal >> external)")

jointPDF_HC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
jointPDF_LC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
sigma_theta = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);

for idx=1:numel(sigma_e)
    modelParams.sigma_i = sqrt(1.5*(sigma_e(idx)^2)); % sigma_i less than sigma_e always
    modelParams.sigma_e = sigma_e(idx);     
    retData = getLLhChoice_MCASANDRE_Contn(modelParams);
    
    jointPDF_HC(idx, :) =  retData.pHC;
    jointPDF_LC(idx, :) =  retData.pLC;
    sigma_theta(idx, :) =  retData.sigma_theta; % sigma_theta is same for all
end


subplot(2, 1, 1)

mean_sigma_theta  = mean(sigma_theta, 2);

hold on

plot(sigma_e', mean_sigma_theta, 'b-', 'LineWidth', 1.5);

xlabel("External noise")
ylabel("Expected [Sigma theta]")

hold off



subplot(2, 1, 2)

hold on

exp_sig_theta_HC = sum(sigma_theta.*jointPDF_HC, 2) ./ sum(jointPDF_HC, 2);
exp_sig_theta_LC = sum(sigma_theta.*jointPDF_LC, 2) ./ sum(jointPDF_LC, 2);

plot(sigma_e', exp_sig_theta_HC, 'b-', 'LineWidth', 1.5, 'DisplayName', 'HC', ...
    DisplayName='HC');
plot(sigma_e', exp_sig_theta_LC, 'r-', 'LineWidth', 1.5, 'DisplayName', 'LC', ...
    DisplayName='LC');

xlabel("External noise")
ylabel("Expected [Sigma theta]")
legend

hold off







sigma_e           = linspace(0.01, 3, 100);

modelParams.sigma_e                   = 0.7;
modelParams.sigma_i                   = 1;
modelParams.varGain                   = 0.5;
modelParams.internalNoiseSamplesCnt   = 200;
modelParams.Cc                        = 1;
modelParams.sigma_m                   = 0.2;

figure
sgtitle("Internal and External noise : Negative covariation")

jointPDF_HC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
jointPDF_LC = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);
sigma_theta = zeros(numel(sigma_e), modelParams.internalNoiseSamplesCnt);

for idx=1:numel(sigma_e)
    sigma_i0 = 0.9;
    k = 0.001;
    modelParams.sigma_i = sqrt( k ./ (sigma_e(idx).^2) + sigma_i0^2 );
    
    modelParams.sigma_e = sigma_e(idx);     
    retData = getLLhChoice_MCASANDRE_Contn(modelParams);
    
    jointPDF_HC(idx, :) =  retData.pHC;
    jointPDF_LC(idx, :) =  retData.pLC;
    sigma_theta(idx, :) =  retData.sigma_theta; % sigma_theta is same for all
end


subplot(2, 1, 1)

mean_sigma_theta  = mean(sigma_theta, 2);

hold on

plot(sigma_e', mean_sigma_theta, 'b-', 'LineWidth', 1.5);

xlabel("External noise")
ylabel("Expected [Sigma theta]")

hold off



subplot(2, 1, 2)

hold on

exp_sig_theta_HC = sum(sigma_theta.*jointPDF_HC, 2) ./ sum(jointPDF_HC, 2);
exp_sig_theta_LC = sum(sigma_theta.*jointPDF_LC, 2) ./ sum(jointPDF_LC, 2);

plot(sigma_e', exp_sig_theta_HC, 'b-', 'LineWidth', 1.5, 'DisplayName', 'HC', ...
    DisplayName='HC');
plot(sigma_e', exp_sig_theta_LC, 'r-', 'LineWidth', 1.5, 'DisplayName', 'LC', ...
    DisplayName='LC');

xlabel("External noise")
ylabel("Expected [Sigma theta]")
legend

hold off


% figure
% 
% for idx=1:numel(sigma_e)
%     subplot(3, 4, idx)
% 
%     x   = sigma_theta(idx, :);
%     pHC = jointPDF_HC(idx, :);
%     pLC = jointPDF_LC(idx, :);
% 
%     hold on
%     plot(x, pHC, DisplayName="HC", LineWidth=1.5);
%     plot(x, pLC, DisplayName="LC", LineWidth=1.5);
%     xlabel("sigma theta")
%     ylabel("Probability")
%     legend
%     hold off
% 
% end


% x = repmat(sigma_e, [modelParams.internalNoiseSamplesCnt 1]); x = x(:); % sigma_e
% y = sigma_theta'; y = y(:); % sigma_theta
% pHC = jointPDF_HC'; pHC = pHC(:);
% pLC = jointPDF_LC'; pLC = pLC(:);
% 
% F = scatteredInterpolant(x, y, pHC, 'natural', 'none');
% [xq, yq] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
% zq = F(xq, yq);
% h1 = imagesc(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100), zq);
% set(h1, 'AlphaData', 0.2);  % Add transparency
% 
% hold on
% 
% F = scatteredInterpolant(x, y, pLC, 'natural', 'none');
% [xq, yq] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
% zq = F(xq, yq);
% h2 = imagesc(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100), zq);
% set(h2, 'AlphaData', 0.2);  % Add transparency
% 
% axis xy;
% xlabel('x'); ylabel('y'); colorbar;
% title('Joint Probability Heatmap');
% 
% hold off


% figure
% % hold on
% scatter3(x, y, pHC, DisplayName="LC");
% scatter3(x, y, pLC, DisplayName="LC");
% xlabel("External noise (contrast level)")
% ylabel("Sigma theta")
% zlabel("pHC")
% % hold off

% imagesc(modelParams.internalNoiseSamples, sigma_e, jointPDF_HC);
% axis xy;  % to make sure the y-axis isn't flipped
% xlabel('Internal noise (\sigma_{int})');
% ylabel('External noise (\sigma_{ext})');
% title('Joint PDF');
% colorbar;

% probConf.HC = zeros(1, numel(sigma_e));
% probConf.LC = zeros(1, numel(sigma_e));
% 
% for idx=1:numel(sigma_e)
% 
%     modelParams.sigma_e                   = sigma_e(idx);
%     modelParams.sigma_i                   = 1;
%     modelParams.varGain                   = 0.5;
%     modelParams.internalNoiseSamplesCnt   = 200;
%     modelParams.Cc                        = 1;
%     modelParams.sigma_m                   = 0.2;
%     
%     probs = getLLhChoice_MCASANDRE_Contn(modelParams);
%     
%     probConf.HC(idx) = probs.pHC;
%     probConf.LC(idx) = probs.pLC;
% end
% 
% figure
% hold on
% plot(sigma_e, probConf.HC, LineWidth=1.5, DisplayName="HC")
% plot(sigma_e, probConf.LC, LineWidth=1.5, DisplayName="LC")
% legend
% xlabel("External noise (contrast level)")
% ylabel("Probability")
% hold off

% % sigma_e = [7, 6, 5, 4, 3, 2, 1]; % In order of increasing contrast level - noise decreases
% sigma_e = [1.9, 1.6, 1.3, 1, 0.7, 0.4, 0.1]; % In order of increasing contrast level - noise decreases
% 
% varGain    = 0.5;
% shapeParam = 1./varGain;
% scaleParam = varGain;
% sigma_i = gamrnd(shapeParam, scaleParam, [1 100]);        % Gain vector for all trials
% 
% figure
% x = 1:numel(sigma_e);
% hold on
% plot(x, sigma_e);
% 
% x1 = repmat(x, [numel(sigma_i) 1]);
% y1 = repmat(sigma_i, [numel(sigma_e) 1])';
% mean_y1 = mean(y1, 1);
% 
% x1 = reshape(x1, [], 1);
% y1 = reshape(y1, [], 1);
% scatter(x1, y1);
% 
% scatter(x, mean_y1, "filled", "o", MarkerEdgeColor='b', MarkerFaceColor='b');
% 
% xlim([0, 7])
% hold off
% 
% 
% 
% sigma_theta = 1; % mu - normal
% sigma_m     = 2; % sigma - normal
% 
% muLogN    = log((sigma_theta.^2) ./ sqrt(sigma_m.^2 + sigma_theta.^2));
% sigmaLogN = sqrt(log((sigma_m.^2)./(sigma_theta.^2) + 1));
% x = linspace(0, 10, 1000);
% pdf = lognpdf(x, muLogN, sigmaLogN);
% samples = lognrnd(muLogN, sigmaLogN, [1 1000]);
% 
% edges = linspace(0, 10, 101);           % 100 bins over [0,10]
% centers = edges(1:end-1) + diff(edges)/2;
% histVals = histcounts(samples, edges, 'Normalization', 'pdf');
% 
% figure
% subplot(1, 2, 1)
% hold on
% bar(centers, histVals, 'hist')          % Manual histogram plot
% plot(x, pdf, 'r-', 'LineWidth', 1.5)    % PDF line
% xlim([0, 10])
% hold off
% 
% % Inverse of log normal distribution
% invSamples = 1./samples;
% 
% sigma_theta = 1; % mu - normal
% sigma_m     = 2; % sigma - normal
% 
% muLogN    = - log((sigma_theta.^2) ./ sqrt(sigma_m.^2 + sigma_theta.^2));
% sigmaLogN = sqrt(log((sigma_m.^2)./(sigma_theta.^2) + 1));
% x = linspace(0, 100, 1000);
% pdf = lognpdf(x, muLogN, sigmaLogN);
% 
% 
% subplot(1, 2, 2)
% hold on
% 
% histogram(invSamples, 'NumBins', 100, 'Normalization', 'pdf');
% plot(x, pdf, 'r-', 'LineWidth', 1.5)
% 
% xlim([0 100])
% hold off


% figure
% hold on
% histogram(samples, NumBins=100, Normalization='probability')
% plot(x, pdf, LineStyle="-", LineWidth=1.5)
% xlim([0, 10])
% hold off