clear all
close all

% Parameters
nTrials = 200;
uniqStimVals = linspace(80, 100, 51);  % Unique set of orientations to show
sigma_d = 2;                          % True stimulus uncertainty
sigma_i = 1;                          % Internal noise
varGain = 0.5;                        % Doubly stochastic sensory noise - Modulated internal noise - Test for different level of variance in gain
Cd = 90;                              % Decision criterion

sigma_m = 0.5;                        % Meta-uncertainty (std of log-normal noise)
Cc = 1.0;                             % Confidence criterion


% Main code starts here 
gainVec = gamrnd(1./varGain, 1./varGain, [numel(uniqStimVals), nTrials]);        % Gain vector for all trials
modulated_sigma_i = sigma_i*gainVec;                                             % Modulated sigma i for each trial
sigma_d_modulated = sqrt(sigma_d^2 + modulated_sigma_i.^2);                      % Modulated stim sigma for each trial

% Psychometric - behavioral performance
sampledStimVal = uniqStimVals' + sigma_d_modulated.*randn(numel(uniqStimVals), nTrials);
choices = zeros(numel(uniqStimVals), nTrials); 
choices(sampledStimVal > Cd) = 1; % CCW choices

psychFn = sum(choices, 2) / size(choices, 2); % Psychometric function

% Plot results
figure
subplot(1, 2, 1)
plot(uniqStimVals, psychFn, LineStyle="--", LineWidth=1.5)
xlabel("Stim value")
ylabel("% CCW")
legend
hold off


subplot(1, 2, 2)
hold on

for idx=1:numel(sigma_m)
    sigma_m_ = sigma_m(idx);

    % Confidence computation
    mu_log = log(sigma_d_modulated.^2 ./ sqrt(sigma_m_^2 + sigma_d_modulated.^2));
    sigma_log = sqrt(log(1 + (sigma_m_^2 ./ sigma_d_modulated.^2)));
    sigma_hat = lognrnd(mu_log, sigma_log);
    
    confidenceVars = abs(sampledStimVal - Cd) ./  sigma_hat;
    confReports = zeros(numel(uniqStimVals), nTrials); 
    confReports(confidenceVars > Cc) = 1; % HC reports
    
    confFn = sum(confReports, 2) / size(confReports, 2); % Confidence function
    
    plot(uniqStimVals, confFn, LineStyle="--", LineWidth=1.5, DisplayName="sigma_m="+ sigma_m_)

end

xlabel("Stim value")
ylabel("% HC")
legend
hold off

% figure
% histogram(sampledStimVal(1, :))

% for i=1:numel(uniqStimVals)
%     stimVal = uniqStimVals(i);
% 
% end

% Plot PDF for each distribution
% x = linspace(-10, 10, 1000);  % X-axis range for plotting
% mu = 0;                       % Assume mean = 0 for all
% 
% figure; hold on;
% for i = 1:length(sigma_stim_modulated)
%     sigma = sigma_stim_modulated(i);
%     y = normpdf(x, mu, sigma);         % Normal distribution with mean 0 and given sigma
%     plot(x, y, 'Color', [0 0 1 0.1]);  % Light blue lines with transparency
% end
% xlabel('x'); ylabel('Probability Density');
% title('Normal Distributions for Each \sigma');
