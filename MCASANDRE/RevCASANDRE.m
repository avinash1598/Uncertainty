clear all
close all

% Parameters
mu_d = linspace(-5, 5, 101);   % Make it an array - Mean of decision variable Vd (e.g., stimulus value)
sigma_d = 1.0;                 % True stimulus uncertainty
sigma_m = [0.5, 10];           % Meta-uncertainty (std of log-normal noise)
Cd = 0.0;                      % Decision criterion
Cc = 0.5;                      % Confidence criterion
n_trials = 500;                % Number of simulated trials

% Preallocate
Vd = zeros(numel(sigma_m), numel(mu_d), n_trials);
estimated_sigma_d = zeros(numel(sigma_m), numel(mu_d), n_trials);
Vc = zeros(numel(sigma_m), numel(mu_d), n_trials);
choice = zeros(numel(sigma_m), numel(mu_d), n_trials);
confidence = zeros(numel(sigma_m), numel(mu_d), n_trials);

% Simulate trials
for jidx=1:numel(sigma_m)
    for sidx=1:length(mu_d)

        stim = mu_d(sidx);
    
        for t = 1:n_trials
            % First stage: simulate noisy estimate of sigma_d (log-normal)
            % Why does it have to be log normal
            mu_log = log(sigma_d^2 / sqrt(sigma_m(jidx)^2 + sigma_d^2));  % log-mean
            sigma_log = sqrt(log(1 + (sigma_m(jidx)^2 / sigma_d^2)));     % log-std
            estimated_sigma_d(jidx, sidx, t) = lognrnd(mu_log, sigma_log);
            
            % estimated_sigma_d(jidx, sidx, t) = sigma_d +
            % sigma_m(jidx)*randn;

            % Second stage: simulate decision variable (normal distribution)
            Vd(jidx, sidx, t) = normrnd(stim, estimated_sigma_d(jidx, sidx, t));
            
            % Make choice
            choice(jidx, sidx, t) = Vd(jidx, sidx, t) > Cd;
        
            % Compute confidence variable (Vc)
            Vc(jidx, sidx, t) = abs(Vd(jidx, sidx, t) - Cd); % / estimated_sigma_d(jidx, sidx, t);
            
            % Generate binary confidence report
            confidence(jidx, sidx, t) = Vc(jidx, sidx, t) > Cc;
        end
    end
end

figure
subplot(2, 4, 1)
hold on

for jidx=1:numel(sigma_m)
    % Psychometric function
    data = squeeze(choice(jidx, :, :));
    psychFn = sum(data, 2) / size(data, 2);

    plot(mu_d, psychFn, LineStyle="--", LineWidth=1.5, DisplayName="sigma_m=" + sigma_m(jidx))
end

xlabel("Stim value")
ylabel("% Choice 1")
legend
hold off

subplot(2, 4, 2)
hold on

for jidx=1:numel(sigma_m)
    % Confidence function
    data = squeeze(confidence(jidx, :, :));
    confFn  = sum(data, 2) / size(data, 2);

    plot(mu_d, confFn, LineStyle="--", LineWidth=1.5, DisplayName="sigma_m=" + sigma_m(jidx))
end

xlabel("Stim value")
ylabel("% HC")
legend
hold off


% subplot(2, 4, 3)
% hold on
% 
% for jidx=1:numel(sigma_m)
%     % Confidence function
%     data = squeeze(Vd(jidx, 51, :));
% 
%     histogram(data, 'NumBins', 20, DisplayName="" + sigma_m(jidx));
% end
% 
% xlabel("Vd")
% ylabel("Count")
% legend
% hold off


% Parameters
mu_d = linspace(-5, 5, 101);   % Make it an array - Mean of decision variable Vd (e.g., stimulus value)
sigma_d = 1.0;                 % True stimulus uncertainty
sigma_m = 0.5;                 % Meta-uncertainty (std of log-normal noise)
Cd = [0.0, 1];                 % Decision criterion
Cc = 0.5;                      % Confidence criterion
n_trials = 500;                % Number of simulated trials

% Preallocate
Vd = zeros(numel(Cd), numel(mu_d), n_trials);
estimated_sigma_d = zeros(numel(Cd), numel(mu_d), n_trials);
Vc = zeros(numel(Cd), numel(mu_d), n_trials);
choice = zeros(numel(Cd), numel(mu_d), n_trials);
confidence = zeros(numel(Cd), numel(mu_d), n_trials);

% Simulate trials
for jidx=1:numel(Cd)
    for sidx=1:length(mu_d)

        stim = mu_d(sidx);
    
        for t = 1:n_trials
            % First stage: simulate noisy estimate of sigma_d (log-normal)
            % Why does it have to be log normal
            mu_log = log(sigma_d^2 / sqrt(sigma_m^2 + sigma_d^2));  % log-mean
            sigma_log = sqrt(log(1 + (sigma_m^2 / sigma_d^2)));     % log-std
            estimated_sigma_d(jidx, sidx, t) = lognrnd(mu_log, sigma_log);
            
            % estimated_sigma_d(jidx, sidx, t) = sigma_d +
            % sigma_m(jidx)*randn;

            % Second stage: simulate decision variable (normal distribution)
            Vd(jidx, sidx, t) = normrnd(stim, estimated_sigma_d(jidx, sidx, t));
            
            % Make choice
            choice(jidx, sidx, t) = Vd(jidx, sidx, t) > Cd(jidx);
        
            % Compute confidence variable (Vc)
            Vc(jidx, sidx, t) = abs(Vd(jidx, sidx, t) - Cd(jidx)); % / estimated_sigma_d(jidx, sidx, t);
            
            % Generate binary confidence report
            confidence(jidx, sidx, t) = Vc(jidx, sidx, t) > Cc;
        end
    end
end

subplot(2, 4, 3)
hold on

for jidx=1:numel(Cd)
    % Psychometric function
    data = squeeze(choice(jidx, :, :));
    psychFn = sum(data, 2) / size(data, 2);

    plot(mu_d, psychFn, LineStyle="--", LineWidth=1.5, DisplayName="Cd=" + Cd(jidx))
end

xlabel("Stim value")
ylabel("% Choice 1")
legend
hold off

subplot(2, 4, 4)
hold on

for jidx=1:numel(Cd)
    % Confidence function
    data = squeeze(confidence(jidx, :, :));
    confFn  = sum(data, 2) / size(data, 2);

    plot(mu_d, confFn, LineStyle="--", LineWidth=1.5, DisplayName="Cd=" + Cd(jidx))
end

xlabel("Stim value")
ylabel("% HC")
legend
hold off



% Parameters
mu_d = linspace(-5, 5, 101);   % Make it an array - Mean of decision variable Vd (e.g., stimulus value)
sigma_d = [1.0, 2.0];          % True stimulus uncertainty
sigma_m = 0.5;                 % Meta-uncertainty (std of log-normal noise)
Cd = 0.0;                      % Decision criterion
Cc = 0.5;                      % Confidence criterion
n_trials = 500;                % Number of simulated trials

% Preallocate
Vd = zeros(numel(sigma_d), numel(mu_d), n_trials);
estimated_sigma_d = zeros(numel(sigma_d), numel(mu_d), n_trials);
Vc = zeros(numel(sigma_d), numel(mu_d), n_trials);
choice = zeros(numel(sigma_d), numel(mu_d), n_trials);
confidence = zeros(numel(sigma_d), numel(mu_d), n_trials);

% Simulate trials
for jidx=1:numel(sigma_d)
    for sidx=1:length(mu_d)

        stim = mu_d(sidx);
    
        for t = 1:n_trials
            % First stage: simulate noisy estimate of sigma_d (log-normal)
            % Why does it have to be log normal
            mu_log = log(sigma_d(jidx)^2 / sqrt(sigma_m^2 + sigma_d(jidx)^2));  % log-mean
            sigma_log = sqrt(log(1 + (sigma_m^2 / sigma_d(jidx)^2)));     % log-std
            estimated_sigma_d(jidx, sidx, t) = lognrnd(mu_log, sigma_log);
            
            % estimated_sigma_d(jidx, sidx, t) = sigma_d +
            % sigma_m(jidx)*randn;

            % Second stage: simulate decision variable (normal distribution)
            Vd(jidx, sidx, t) = normrnd(stim, estimated_sigma_d(jidx, sidx, t));
            
            % Make choice
            choice(jidx, sidx, t) = Vd(jidx, sidx, t) > Cd;
        
            % Compute confidence variable (Vc)
            Vc(jidx, sidx, t) = abs(Vd(jidx, sidx, t) - Cd); % / estimated_sigma_d(jidx, sidx, t);
            
            % Generate binary confidence report
            confidence(jidx, sidx, t) = Vc(jidx, sidx, t) > Cc;
        end
    end
end

subplot(2, 4, 5)
hold on

for jidx=1:numel(sigma_d)
    % Psychometric function
    data = squeeze(choice(jidx, :, :));
    psychFn = sum(data, 2) / size(data, 2);

    plot(mu_d, psychFn, LineStyle="--", LineWidth=1.5, DisplayName="sigma_d=" + sigma_d(jidx))
end

xlabel("Stim value")
ylabel("% Choice 1")
legend
hold off

subplot(2, 4, 6)
hold on

for jidx=1:numel(sigma_d)
    % Confidence function
    data = squeeze(confidence(jidx, :, :));
    confFn  = sum(data, 2) / size(data, 2);

    plot(mu_d, confFn, LineStyle="--", LineWidth=1.5, DisplayName="sigma_d=" + sigma_d(jidx))
end

xlabel("Stim value")
ylabel("% HC")
legend
hold off



% Parameters
mu_d = linspace(-5, 5, 101);   % Make it an array - Mean of decision variable Vd (e.g., stimulus value)
sigma_d = 1.0;                 % True stimulus uncertainty
sigma_m = 0.5;                 % Meta-uncertainty (std of log-normal noise)
Cd = 0.0;                      % Decision criterion
Cc = [0.5, 1.5];               % Confidence criterion
n_trials = 500;                % Number of simulated trials

% Preallocate
Vd = zeros(numel(Cc), numel(mu_d), n_trials);
estimated_sigma_d = zeros(numel(Cc), numel(mu_d), n_trials);
Vc = zeros(numel(Cc), numel(mu_d), n_trials);
choice = zeros(numel(Cc), numel(mu_d), n_trials);
confidence = zeros(numel(Cc), numel(mu_d), n_trials);

% Simulate trials
for jidx=1:numel(Cc)
    for sidx=1:length(mu_d)

        stim = mu_d(sidx);
    
        for t = 1:n_trials
            % First stage: simulate noisy estimate of sigma_d (log-normal)
            % Why does it have to be log normal
            mu_log = log(sigma_d^2 / sqrt(sigma_m^2 + sigma_d^2));  % log-mean
            sigma_log = sqrt(log(1 + (sigma_m^2 / sigma_d^2)));     % log-std
            estimated_sigma_d(jidx, sidx, t) = lognrnd(mu_log, sigma_log);
            
            % estimated_sigma_d(jidx, sidx, t) = sigma_d +
            % sigma_m(jidx)*randn;

            % Second stage: simulate decision variable (normal distribution)
            Vd(jidx, sidx, t) = normrnd(stim, estimated_sigma_d(jidx, sidx, t));
            
            % Make choice
            choice(jidx, sidx, t) = Vd(jidx, sidx, t) > Cd;
        
            % Compute confidence variable (Vc)
            Vc(jidx, sidx, t) = abs(Vd(jidx, sidx, t) - Cd); % / estimated_sigma_d(jidx, sidx, t);
            
            % Generate binary confidence report
            confidence(jidx, sidx, t) = Vc(jidx, sidx, t) > Cc(jidx);
        end
    end
end

subplot(2, 4, 7)
hold on

for jidx=1:numel(Cc)
    % Psychometric function
    data = squeeze(choice(jidx, :, :));
    psychFn = sum(data, 2) / size(data, 2);

    plot(mu_d, psychFn, LineStyle="--", LineWidth=1.5, DisplayName="Cc=" + Cc(jidx))
end

xlabel("Stim value")
ylabel("% Choice 1")
legend
hold off

subplot(2, 4, 8)
hold on

for jidx=1:numel(Cc)
    % Confidence function
    data = squeeze(confidence(jidx, :, :));
    confFn  = sum(data, 2) / size(data, 2);

    plot(mu_d, confFn, LineStyle="--", LineWidth=1.5, DisplayName="Cc=" + Cc(jidx))
end

xlabel("Stim value")
ylabel("% HC")
legend
hold off


