clear all
close all

% This is probably not right!

niterations = 100;
ntrials_per_ori_list = linspace(100, 2000, 20);

resultsBootstrapTest = zeros(numel(ntrials_per_ori_list), 1);
resultsBootstrapTrendTest = zeros(numel(ntrials_per_ori_list), 1);
resultsBootstrapTest_ValidItr = zeros(numel(ntrials_per_ori_list), 1);
resultsBootstrapTrendTest_ValidItr = zeros(numel(ntrials_per_ori_list), 1);

for tidx=1:numel(ntrials_per_ori_list)

    tempBootstrap = 0;
    tempBootstrapTrend = 0;
    totalValidCnt = 0;

    for itr=1:niterations
        orientations     = linspace(0, 180, 2); % linspace(0, 180, 18);
        ntrials_per_ori  = ntrials_per_ori_list(tidx);
        sigma_e_levels   = linspace(0.01, 2, 4);  % Varying external noise or contrast
        % sigma_i          = 1;
        varGain          = 0.5;
        sigma_m          = 0.2;
        Cc               = 0.5;
        
        % Preallocate arrays
        n_theta                  = numel(orientations);
        n_levels                 = numel(sigma_e_levels);
        
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
                sigma_i = sqrt(9*(sigma_e^2));
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
        
        % IMP: All analysis done for single orientation
        
        % Plot the performance curves
        resp_err_all = abs(theta_resp_all - theta_true_all);
        % resp_err_all_reshaped = reshape(resp_err_all, numel(sigma_e_levels), []);
        % confidence_report_all_reshaped = reshape(confidence_report_all, numel(sigma_e_levels), []);
        
        resp_err_all_reshaped = squeeze(resp_err_all(:,1,:)); % Data for one angle only
        confidence_report_all_reshaped = squeeze(confidence_report_all(:,1,:));
        
        % Bootstrap test
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
        x3 = resp_HC(numel(sigma_e_levels), :); valid_idx = ~isnan(x3); x3 = x3(valid_idx);
        x4 = resp_LC(numel(sigma_e_levels), :); valid_idx = ~isnan(x4); x4 = x4(valid_idx);
    
        if (numel(x1) > 0) && (numel(x2) > 0) && (numel(x3) > 0) && (numel(x4) > 0)
            % Bootstrap test
            result = doBootstrapTest(x1, x2, x3, x4, true);
            tempBootstrap = tempBootstrap + (result.pval < 0.05);

            totalValidCnt = totalValidCnt + 1;
        end

        % Bootstrap trend test
        result = doBootstrapTrendTest(sigma_e_levels, resp_HC, resp_LC, true);
        tempBootstrapTrend = tempBootstrapTrend + (result.pval < 0.05);    
    end
    
    resultsBootstrapTest(tidx) = tempBootstrap / niterations;
    resultsBootstrapTrendTest(tidx) = tempBootstrapTrend / niterations;
    
    resultsBootstrapTest_ValidItr(tidx) = tempBootstrap / totalValidCnt;
    resultsBootstrapTrendTest_ValidItr(tidx) = tempBootstrapTrend / totalValidCnt;
    
    if niterations ~= totalValidCnt
        fprintf("Mismatch \n")
    end
end


figure 

subplot(2, 2, 1)
plot(ntrials_per_ori_list, resultsBootstrapTest)
xlabel("Trials per orientation")
ylabel("% significant")
title("Bootstrap test")

subplot(2, 2, 3)
plot(ntrials_per_ori_list, resultsBootstrapTest_ValidItr)
xlabel("Trials per orientation")
ylabel("% significant")
title("Bootstrap test (Iterations having both HC and LC reports)")

subplot(2, 2, 2)
plot(ntrials_per_ori_list, resultsBootstrapTrendTest)
xlabel("Trials per orientation")
ylabel("% significant")
title("Bootstrap trend test")

subplot(2, 2, 4)
plot(ntrials_per_ori_list, resultsBootstrapTrendTest_ValidItr)
xlabel("Trials per orientation")
ylabel("% significant")
title("Bootstrap trend test (Iterations having both HC and LC reports")


