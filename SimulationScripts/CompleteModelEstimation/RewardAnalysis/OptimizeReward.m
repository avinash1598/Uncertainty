clear all
close all

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
    end
end


%% Fit model
metaData.theta_true_all  = theta_true_all(:);
metaData.theta_resp_all  = theta_resp_all(:);
metaData.VcEstimates_all = VcEstimates(:);

p_maxTolerableError = 25;
p_y1                = 10;
p_c1                = 0.2;
p_m1_penalty        = 1;
p_penalty_HC        = 0.2;

params = [p_maxTolerableError p_y1 p_c1 p_m1_penalty p_penalty_HC];

% Objective function
objFun = @(x) minimizeError(x, metaData);

% Bounds (ga requires finite bounds!)
lb = zeros(size(params));     % same as before
ub = []; % example finite upper bounds

% Optimization options for fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...          % or 'interior-point', 'trust-region-reflective', etc.
    'MaxIterations', 1000, ...
    'OptimalityTolerance', 1e-20, ...
    'StepTolerance', 1e-20);

% Initial guess (required for fmincon)
x0 = params;   % start in the middle of bounds, for example

% Run fmincon
[optimalValues, fval, exitflag, output] = fmincon(objFun, x0, ...
    [], [], [], [], lb, ub, [], options);

disp('Optimal parameters:');
disp(optimalValues);
disp('Final objective value:');
disp(fval);


%%
function loss = minimizeError(params, data)

theta_true_all = data.theta_true_all;
theta_resp_all = data.theta_resp_all;
VcEstimates_all = data.VcEstimates_all;

% Ccs                        = linspace(0, 0.5, 200);
% HC_tot_reward_List         = zeros(size(Ccs));
% LC_tot_reward_List         = zeros(size(Ccs));
% tot_reward_HC_LC           = zeros(size(Ccs));
% total_reward_all_LC        = zeros(size(Ccs));
% total_reward_all_HC        = zeros(size(Ccs));
% total_reward_random_guess  = zeros(size(Ccs));

% for i=1:numel(Ccs)
    Cc_ = 0.06;
    confidence_report_all_  = VcEstimates_all > Cc_;
    rewardVals = calcReward(theta_true_all, theta_resp_all, confidence_report_all_, params);

    totalReward_HC = sum(rewardVals(confidence_report_all_ == 1));
    totalReward_LC = sum(rewardVals(confidence_report_all_ == 0));

    tot_reward_HC_LC = totalReward_HC + totalReward_LC;
    HC_tot_reward_List = totalReward_HC;
    LC_tot_reward_List = totalReward_LC;

    % All HC
    conf_report = 1 + zeros(size(confidence_report_all_));
    rewardVals = calcReward(theta_true_all, theta_resp_all, conf_report, params);
    total_reward_all_HC = sum(rewardVals);

    % All LC
    conf_report = zeros(size(confidence_report_all_));
    rewardVals = calcReward(theta_true_all, theta_resp_all, conf_report, params);
    total_reward_all_LC = sum(rewardVals);
    
    % Random guesses
    randVals = rand(size(confidence_report_all_));
    conf_report = randVals > 0.5;
    rewardVals = calcReward(theta_true_all, theta_resp_all, conf_report, params);
    total_reward_random_guess = sum(rewardVals);

% end

% Minimize difference between total_reward_random_guess and zero, HC or LC only
% Maximize diff between total_reward_random_guess and optimal total reward
loss = ( total_reward_random_guess - 0 ).^2 + ...
        ( total_reward_random_guess - total_reward_all_LC ).^2 + ...
        ( total_reward_random_guess - total_reward_all_HC ).^2 ;

%- ...
%        (total_reward_random_guess - tot_reward_HC_LC).^2;

loss = sum(loss);

end

% Calculate reward
function reward = calcReward(trueOri, reportedOri, confReport, params)

maxTolerableError = params(1);
y1                = params(2);
c1                = params(3);
m1_penalty        = params(4);
penalty_HC        = params(5);

% Note: reported orientation is already pi-periodic, true orientation as
% well
% maxTolerableError = maxTolerableError; % In degrees
absPerceptualError = abs(trueOri - reportedOri);
absPerceptualError = min(absPerceptualError, 180 - absPerceptualError);
reward = zeros(size(absPerceptualError));

% y1 = 10;   % variable
x1 = 9;    % fixed
% c1 = 0.2;  % variable
m1 = - c1/y1;
m2 =  ( m1*x1 + c1 ) / (x1 - maxTolerableError); % c2 / maxTolerableError; %m1 - (c1 - c2)/x1;
c2 = - m2*maxTolerableError;

% High confidence
gHC = m1 * absPerceptualError + c1;
reward(:) = gHC;
idx = (absPerceptualError <= y1) & (confReport == 1);
reward(idx) = gHC(idx);
idx = (absPerceptualError > y1) & (confReport == 1);
reward(idx) = m1_penalty*gHC(idx);

% Low confidence
gLC = m2 * absPerceptualError + c2;
reward(confReport == 0) = gLC(confReport == 0);

% HC error
idx = ( absPerceptualError > maxTolerableError ) & (confReport == 1);
reward(idx) = - penalty_HC; %-0.4; 

% LC error
idx = ( absPerceptualError > maxTolerableError ) & (confReport == 0);
reward(idx) = 0;

end
