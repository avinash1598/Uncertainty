function data = GenerateIndData(modelParams, n_uncertainty_levels, modelType, orientations, ntrials_per_ori)
    
    assert(modelType == "reduced" || modelType=="full" || modelType == "full_jumbo")
    
    if modelType == "reduced"
        assert(numel(modelParams) == 11);
    elseif modelType == "full"
        assert(numel(modelParams) == 13);
    elseif modelType == "full_jumbo"
        assert(numel(modelParams) == 18);
    end

    %orientations     = 0:15:175; %linspace(0, 180, 18); %0:10:180; % linspace(0, 180, 18);
    %ntrials_per_ori  = 24;
    
    sigma_s        = modelParams(1:n_uncertainty_levels);
    shape          = modelParams(n_uncertainty_levels + 1);
    scale          = modelParams(n_uncertainty_levels + 2);
    sigma_meta     = modelParams(n_uncertainty_levels + 3);
    Cc             = modelParams(n_uncertainty_levels + 4);
    guessRate      = modelParams(n_uncertainty_levels + 5);
    
    % Preallocate arrays
    n_theta                  = numel(orientations);
    theta_true_all           = zeros(n_uncertainty_levels, n_theta, ntrials_per_ori);
    theta_resp_all           = zeros(n_uncertainty_levels, n_theta, ntrials_per_ori); % Recorded theta based on user response
    confidence_report_all    = zeros(n_uncertainty_levels, n_theta, ntrials_per_ori);
    
    if modelType == "reduced"
        bias            = zeros(1, numel(orientations));
        sigma_s_stim    = repmat(sigma_s', [1 n_theta]);
        bias(:)         = 0;
    elseif modelType == "full"
        sigma_ori_scale = modelParams(n_uncertainty_levels + 6);
        biasAmp         = modelParams(n_uncertainty_levels + 7);
        sigma_s_stim    = sigma_s' + (sigma_ori_scale.*sigma_s)'*(abs(sind(orientations - 90)));
        bias            = biasAmp*sind(2*orientations); 
    elseif modelType == "full_jumbo"
        sigma_a          = modelParams(n_uncertainty_levels + 6:2*n_uncertainty_levels + 5);
        biasAmp          = modelParams(2*n_uncertainty_levels + 6);
        sigma_s_stim     = sigma_s' + (sigma_a)'*(abs(sind(orientations - 90)));
        bias             = biasAmp*sind(2*orientations); 
    end
    
    for l=1:n_uncertainty_levels
        for i = 1:n_theta
            theta_true = orientations(i);   % True orientation
            trials = ntrials_per_ori;
            
            % Step1: Using estimate of this uncertainty, the subject estimates
            % the orientation
            % Internal estimate (Gaussian noise) - Note: this is not wraped
            % In actual behavioral data this will be wraped
            % theta_est = theta_true + sigma_p * randn(trials, 1);
            % Find doubly stochastic theta
            gain = gamrnd(shape, scale, [trials 1]);
            sigma_si_modulated = gain;
            sigma_m_stim = sqrt(sigma_s_stim(l, i).^2 + sigma_si_modulated);
            mean_m_stim = theta_true + bias(i);
            
            % TODO: take into account bias?
            % Wrap the angle at the plotting stage. Note: warapping should be
            % performed according to the true angle.
            theta_est = mean_m_stim + sigma_m_stim .* randn(trials, 1);
            theta_est = mod(theta_est, 180); 
            
            % Simulat guess rate
            guess_tl_idx = randi([1 trials], floor( trials*guessRate ), 1);
            guessOris = 180*rand(numel(guess_tl_idx), 1);
            theta_est(guess_tl_idx) = guessOris;
            
            assert(numel(sigma_m_stim) == trials);
            
            % Step1: Subject first gets an estimate of its uncertainty
            % Subject’s estimate of their uncertainty (meta-uncertainty)
            mu_log = log(sigma_m_stim.^2 ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2));
            sigma_log = sqrt(log(1 + (sigma_meta.^2 ./ sigma_m_stim.^2)));
            sigma_hat = lognrnd(mu_log, sigma_log, trials, 1);
            
            % Confidence variable
            Vc = 1 ./ sigma_hat;
            
            % Store
            theta_true_all(l, i, :)         = theta_true;
            theta_resp_all(l, i, :)         = theta_est;
            confidence_report_all(l, i, :)  = Vc > Cc;
        end
    end

    % Prepare data structures 
    resp_err_all = (theta_resp_all - theta_true_all);
    resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error
    
    % Save model data
    data.theta_true_all         = theta_true_all;
    data.theta_resp_all         = theta_resp_all;
    data.resp_err_all           = resp_err_all;
    data.confidence_report_all  = confidence_report_all;
    data.stdByOri               = squeeze( std(resp_err_all, 0, 3) );
    data.madByOri               = squeeze( mad(resp_err_all, 1, 3) );
    data.orientations           = orientations';
    
    T = table( ...
        data.theta_true_all(:), ...
        data.resp_err_all(:), ...
        'VariableNames', {'theta_true', 'resp_err'});
    G = groupsummary(T, {'theta_true'}, {'mean'}, 'resp_err');
    data.bias = G.mean_resp_err;

    data.params.sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 )';
    %data.params.b                     = sigma_s;
    %data.params.a                     = sigma_ori_scale*sigma_s;
    %data.params.biasAmp               = biasAmp;
    data.params.shape                 = shape;
    data.params.scale                 = scale;
    data.params.sigma_meta            = sigma_meta;
    data.params.Cc                    = Cc;
    data.params.guessRate             = guessRate;

end