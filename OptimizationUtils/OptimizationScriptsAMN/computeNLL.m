% Loss function for optimization
function nll = computeNLL(params, metaData)

% Metadata
nLevels              = metaData.n_levels;
trlErrors            = metaData.trlErrors;
trlConfReports       = metaData.trlConfReports;
trlUncertaintyLevels = metaData.trlUncertaintyLevels;

% Params
param_sigma_ext          = params(1:nLevels);
param_ext_noise_scale    = params(nLevels + 1:2*nLevels);
param_scale              = params(2*nLevels + 1);
param_sigma_meta         = params(2*nLevels + 2);
param_Cc                 = params(2*nLevels + 3);
param_guessrate          = params(2*nLevels + 4);

trial_probs      = zeros(size(trlErrors));
trial_probs_HC   = zeros(size(trlErrors));   % conditional prob of error given conf report - HC
trial_probs_LC   = zeros(size(trlErrors));   % conditional prob of error given conf report - LC
trial_probs_Conf = zeros(size(trlErrors)); % Probability of a trial being high or low confidence

for i=1:nLevels
    modelParams.sigma_ext           = param_sigma_ext(i);
    modelParams.ext_noise_scale     = param_ext_noise_scale(i);
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = param_guessrate;
    %modelParams.oriErrBinWidth      = 0.1;
    
    retData = getPDFs(modelParams, true); % originally set to true
    
    % PDF
    idx = (trlUncertaintyLevels == i);
    trial_probs(idx) = interp1( ...
        retData.rvOriErrs, ...
        retData.analyticalPDF(:), ...
        trlErrors(idx), ...
        'linear');
    
    % PDF HC
    idxHC = (trlConfReports == 1) & idx;
    trial_probs_HC(idxHC) = interp1( ...
        retData.rvOriErrs, ...
        retData.analyticalPDF_HC(:), ...
        trlErrors(idxHC), ...
        'linear');
    trial_probs_Conf(idxHC) = retData.pHC;
    
    % PDF LC
    idxLC = (trlConfReports == 0) & idx;
    trial_probs_LC(idxLC) = interp1( ...
        retData.rvOriErrs, ...
        retData.analyticalPDF_LC(:), ...
        trlErrors(idxLC), ...
        'linear');
    trial_probs_Conf(idxLC) = retData.pLC; % Mutually exclusive from LC
    
end

% NLL loss
% P(err|trial) = P(err|trial, confreport)*P(confreport|trial)
ll_HC = log( trial_probs_HC .* trial_probs_Conf + eps); % P(ERR|HC)*P(HC)
ll_LC = log( trial_probs_LC .* trial_probs_Conf + eps); % P(ERR|LC)*P(LC)
% ll    = log( trial_probs + eps ); 

% ll_HC and ll_LC intersection should be zero. Or is it?

nll = ( ll_HC + ll_LC); %ll +  % Summation of log => correct thing to do
nll = - sum(nll(:));

end
