
% clear all
% close all

% Parameters
nTrials = 500;
uniqStimVals = linspace(80, 100, 51); % Unique set of orientations to show
sigma_d = [2, 10];                     % True stimulus uncertainty
sigma_i = 1;                          % Internal noise
varGain = 0.5;                        % Doubly stochastic sensory noise - Modulated internal noise
Cd = 90;                              % Decision criterion

sigma_m = 0.5;                        % Meta-uncertainty (std of log-normal noise)
Cc = 1.0;                             % Confidence criterion


% Plot results
figure
hold on

for idx=1:numel(sigma_d)
    sigma_d_ = sigma_d(idx);

    % Main code starts here 
    gainVec = gamrnd(1./varGain, 2, [numel(uniqStimVals), nTrials]);        % Gain vector for all trials
    modulated_sigma_i = sigma_i*gainVec;                           % Modulated sigma i for each trial
    sigma_d_modulated = sqrt(sigma_d_^2 + modulated_sigma_i.^2); % Modulated stim sigma for each trial
    
    % Psychometric - behavioral performance
    sampledStimVal = uniqStimVals' + sigma_d_modulated.*randn(numel(uniqStimVals), nTrials);
    choices = zeros(numel(uniqStimVals), nTrials); 
    choices(sampledStimVal > Cd) = 1; % CCW choices
    
    psychFn = sum(choices, 2) / size(choices, 2); % Psychometric function

    % Confidence computation
    mu_log = log(sigma_d_modulated.^2 ./ sqrt(sigma_m^2 + sigma_d_modulated.^2));
    sigma_log = sqrt(log(1 + (sigma_m^2 ./ sigma_d_modulated.^2)));
    sigma_hat = lognrnd(mu_log, sigma_log);
    
    confidenceVarsVals = (sampledStimVal - Cd) ./  sigma_hat;
    confidenceVars = abs(sampledStimVal - Cd) ./  sigma_hat;
    confReports = zeros(numel(uniqStimVals), nTrials); 
    confReports(confidenceVars > Cc) = 1; % HC reports
    
    confFn = sum(confReports, 2) / size(confReports, 2); % Confidence function


    subplot(2, 2, 1); hold on;
    plot(uniqStimVals, psychFn, LineStyle="-", LineWidth=1.5, DisplayName="sigma_d="+ sigma_d_)
    xlabel("Stim value")
    ylabel("% CCW")
    legend

    subplot(2, 2, 2); hold on;
    plot(uniqStimVals, confFn, LineStyle="-", LineWidth=1.5, DisplayName="sigma_d="+ sigma_d_)
    xlabel("Stim value")
    ylabel("% HC")
    legend

    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 51;
    data = sampledStimVal(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 3); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="sigma_d="+ sigma_d_);
    xline(Cd, LineStyle="--", HandleVisibility='off');
    xlabel('Stim value');
    title('PDF');
    legend

    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 51;
    data = confidenceVarsVals(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 4); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="sigma_d="+ sigma_d_);
    xline(Cc, LineStyle="--", HandleVisibility='off');
    xline(-Cc, LineStyle="--", HandleVisibility='off');
    xlabel('Confidence Var');
    title('PDF');
    legend
end

hold off

sgtitle('Change in sigma_d');



% Parameters
nTrials = 200;
uniqStimVals = linspace(80, 100, 51); % Unique set of orientations to show
sigma_d = 2;                          % True stimulus uncertainty
sigma_i = 1;                          % Internal noise
varGain = 0.5;                        % Doubly stochastic sensory noise - Modulated internal noise
Cd = [90, 93];                        % Decision criterion

sigma_m = 0.5;                        % Meta-uncertainty (std of log-normal noise)
Cc = 1.0;                             % Confidence criterion


% Plot results
figure
hold on

for idx=1:numel(Cd)
    Cd_ = Cd(idx);

    % Main code starts here 
    gainVec = gamrnd(1./varGain, 2, [numel(uniqStimVals), nTrials]);        % Gain vector for all trials
    modulated_sigma_i = sigma_i*gainVec;                           % Modulated sigma i for each trial
    sigma_d_modulated = sqrt(sigma_d^2 + modulated_sigma_i.^2); % Modulated stim sigma for each trial
    
    % Psychometric - behavioral performance
    sampledStimVal = uniqStimVals' + sigma_d_modulated.*randn(numel(uniqStimVals), nTrials);
    choices = zeros(numel(uniqStimVals), nTrials); 
    choices(sampledStimVal > Cd_) = 1; % CCW choices
    
    psychFn = sum(choices, 2) / size(choices, 2); % Psychometric function

    % Confidence computation
    mu_log = log(sigma_d_modulated.^2 ./ sqrt(sigma_m^2 + sigma_d_modulated.^2));
    sigma_log = sqrt(log(1 + (sigma_m^2 ./ sigma_d_modulated.^2)));
    sigma_hat = lognrnd(mu_log, sigma_log);
    
    confidenceVarsVals = (sampledStimVal - Cd_) ./  sigma_hat;
    confidenceVars = abs(sampledStimVal - Cd_) ./  sigma_hat;
    confReports = zeros(numel(uniqStimVals), nTrials); 
    confReports(confidenceVars > Cc) = 1; % HC reports
    
    confFn = sum(confReports, 2) / size(confReports, 2); % Confidence function


    subplot(2, 2, 1); hold on;
    plot(uniqStimVals, psychFn, LineStyle="-", LineWidth=1.5, DisplayName="Cd="+ Cd_)
    xlabel("Stim value")
    ylabel("% CCW")
    legend

    subplot(2, 2, 2); hold on;
    plot(uniqStimVals, confFn, LineStyle="-", LineWidth=1.5, DisplayName="Cd="+ Cd_)
    xlabel("Stim value")
    ylabel("% HC")
    legend
    
    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 51;
    data = sampledStimVal(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 3); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="Cd="+ Cd_);
    xline(Cd, LineStyle="--", HandleVisibility='off');
    xlabel('Stim value');
    title('PDF');
    legend

    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 51;
    data = confidenceVarsVals(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 4); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="Cd="+ Cd_);
    xline(Cc, LineStyle="--", HandleVisibility='off');
    xline(-Cc, LineStyle="--", HandleVisibility='off');
    xlabel('Confidence Var');
    title('PDF');
    legend
end

hold off

sgtitle('Change in Cd');





% Parameters
nTrials = 200;
uniqStimVals = linspace(45, 135, 101); % Unique set of orientations to show
sigma_d = 2;                          % True stimulus uncertainty
sigma_i = 1;                          % Internal noise
varGain = 0.5;                        % Doubly stochastic sensory noise - Modulated internal noise
Cd = 90;                              % Decision criterion

sigma_m = [0.5, 10];                   % Meta-uncertainty (std of log-normal noise)
Cc = 1.0;                             % Confidence criterion


% Plot results
figure
hold on

for idx=1:numel(sigma_m)
    sigma_m_ = sigma_m(idx);

    % Main code starts here 
    gainVec = gamrnd(1./varGain, 2, [numel(uniqStimVals), nTrials]);        % Gain vector for all trials
    modulated_sigma_i = sigma_i*gainVec;                           % Modulated sigma i for each trial
    sigma_d_modulated = sqrt(sigma_d^2 + modulated_sigma_i.^2); % Modulated stim sigma for each trial
    
    % Psychometric - behavioral performance
    sampledStimVal = uniqStimVals' + sigma_d_modulated.*randn(numel(uniqStimVals), nTrials);
    choices = zeros(numel(uniqStimVals), nTrials); 
    choices(sampledStimVal > Cd) = 1; % CCW choices
    
    psychFn = sum(choices, 2) / size(choices, 2); % Psychometric function

    % Confidence computation
    mu_log = log(sigma_d_modulated.^2 ./ sqrt(sigma_m_^2 + sigma_d_modulated.^2));
    sigma_log = sqrt(log(1 + (sigma_m_^2 ./ sigma_d_modulated.^2)));
    sigma_hat = lognrnd(mu_log, sigma_log);
    
    confidenceVarsVals = (sampledStimVal - Cd) ./  sigma_hat;
    confidenceVars = abs(sampledStimVal - Cd) ./  sigma_hat;
    confReports = zeros(numel(uniqStimVals), nTrials); 
    confReports(confidenceVars > Cc) = 1; % HC reports
    
    confFn = sum(confReports, 2) / size(confReports, 2); % Confidence function


    subplot(2, 2, 1); hold on;
    plot(uniqStimVals, psychFn, LineStyle="-", LineWidth=1.5, DisplayName="sigma_m="+ sigma_m_)
    xlabel("Stim value")
    ylabel("% CCW")
    legend

    subplot(2, 2, 2); hold on;
    plot(uniqStimVals, confFn, LineStyle="-", LineWidth=1.5, DisplayName="sigma_m="+ sigma_m_)
    xlabel("Stim value")
    ylabel("% HC")
    legend
    
    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 71;
    data = sampledStimVal(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 3); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="sigma_m="+ sigma_m_);
    xline(Cd, LineStyle="--", HandleVisibility='off');
    xlabel('Stim value');
    title('PDF');
    legend
    
    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 71;
    data = confidenceVarsVals(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 4); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="sigma_m="+ sigma_m_);
    xline(Cc, LineStyle="--", HandleVisibility='off');
    xline(-Cc, LineStyle="--", HandleVisibility='off');
    xlabel('Confidence Var');
    title('PDF');
    legend
end

hold off
sgtitle('Change in sigma_m');




% Parameters
nTrials = 200;
uniqStimVals = linspace(80, 100, 51); % Unique set of orientations to show
sigma_d = 2;                          % True stimulus uncertainty
sigma_i = 1;                          % Internal noise
varGain = 0.5;                        % Doubly stochastic sensory noise - Modulated internal noise
Cd = 90;                              % Decision criterion

sigma_m = 0.5;                        % Meta-uncertainty (std of log-normal noise)
Cc = [1.0, 2.0];                      % Confidence criterion


% Plot results
figure
hold on

for idx=1:numel(Cc)
    Cc_ = Cc(idx);

    % Main code starts here 
    gainVec = gamrnd(1./varGain, 2, [numel(uniqStimVals), nTrials]);        % Gain vector for all trials
    modulated_sigma_i = sigma_i*gainVec;                           % Modulated sigma i for each trial
    sigma_d_modulated = sqrt(sigma_d^2 + modulated_sigma_i.^2); % Modulated stim sigma for each trial
    
    % Psychometric - behavioral performance
    sampledStimVal = uniqStimVals' + sigma_d_modulated.*randn(numel(uniqStimVals), nTrials);
    choices = zeros(numel(uniqStimVals), nTrials); 
    choices(sampledStimVal > Cd) = 1; % CCW choices
    
    psychFn = sum(choices, 2) / size(choices, 2); % Psychometric function

    % Confidence computation
    mu_log = log(sigma_d_modulated.^2 ./ sqrt(sigma_m^2 + sigma_d_modulated.^2));
    sigma_log = sqrt(log(1 + (sigma_m^2 ./ sigma_d_modulated.^2)));
    sigma_hat = lognrnd(mu_log, sigma_log);
    
    confidenceVarsVals = (sampledStimVal - Cd) ./  sigma_hat;
    confidenceVars = abs(sampledStimVal - Cd) ./  sigma_hat;
    confReports = zeros(numel(uniqStimVals), nTrials); 
    confReports(confidenceVars > Cc_) = 1; % HC reports
    
    confFn = sum(confReports, 2) / size(confReports, 2); % Confidence function


    subplot(2, 2, 1); hold on;
    plot(uniqStimVals, psychFn, LineStyle="-", LineWidth=1.5, DisplayName="Cc="+ Cc_)
    xlabel("Stim value")
    ylabel("% CCW")
    legend

    subplot(2, 2, 2); hold on;
    plot(uniqStimVals, confFn, LineStyle="-", LineWidth=1.5, DisplayName="Cc="+ Cc_)
    xlabel("Stim value")
    ylabel("% HC")
    legend
    
    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 51;
    data = sampledStimVal(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 3); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="Cc="+ Cc_);
    xline(Cd, LineStyle="--", HandleVisibility='off');
    xlabel('Stim value');
    title('PDF');
    legend

    % Pick one of the stim and show sampled value and confidence
    % disribution
    stimIdx = 51;
    data = confidenceVarsVals(stimIdx, :);
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;
    
    subplot(2, 2, 4); hold on;
    plot(binCenters, counts, LineStyle="-", LineWidth=1.5, DisplayName="Cc="+ Cc_);
    xline(Cc, LineStyle="--", HandleVisibility='off');
    xline(-Cc, LineStyle="--", HandleVisibility='off');
    xlabel('Confidence Var');
    xlim([-10, 10])
    title('PDF');
    legend
end

hold off
sgtitle('Change in Cc');
