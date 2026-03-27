% Loss function for optimization
function nll = computeNLL(params, metaData)

% Metadata
nLevels              = metaData.n_levels;
trlErrors            = metaData.trlErrors;
trlConfReports       = metaData.trlConfReports;
trlUncertaintyLevels = metaData.trlUncertaintyLevels;
trlStimOris          = metaData.trlStimOris;
stimulusEnergy       = metaData.stimulusEnergy;
sigma_m_shape1       = metaData.sigma_m_shape1;
sigma_m_shape2       = metaData.sigma_m_shape2;
biasShape            = metaData.biasShape;   
orientations         = metaData.orientations;

% Params
param_sigma_ext          = params(1:nLevels);
param_b                  = params(nLevels+1:2*nLevels);
param_scale              = params(2*nLevels + 1);
param_sigma_meta         = params(2*nLevels + 2);
param_Cc                 = params(2*nLevels + 3);
param_guessrate          = params(2*nLevels + 4);
param_ori_scale          = params(2*nLevels + 5);
param_bias               = params(2*nLevels + 6);

trial_probs      = zeros(size(trlErrors));
trial_probs_HC   = zeros(size(trlErrors));   % conditional prob of error given conf report - HC
trial_probs_LC   = zeros(size(trlErrors));   % conditional prob of error given conf report - LC
trial_probs_Conf = zeros(size(trlErrors)); % Probability of a trial being high or low confidence

for i=1:nLevels
    modelParams.sigma_ext           = param_sigma_ext(i);
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = param_guessrate; % constant guess rate
    modelParams.stimulusEnergy      = stimulusEnergy(i);
    %modelParams.oriErrBinWidth      = 0.1;
    
    modelParams.b              = param_b(i);
    modelParams.a              = param_ori_scale*param_b(i);
    modelParams.biasAmp        = param_bias;
    modelParams.sigma_m_shape1 = sigma_m_shape1;
    modelParams.sigma_m_shape2 = sigma_m_shape2;
    modelParams.biasShape      = biasShape;
    
    retData = getPDFs_v4(modelParams, orientations, true); % originally set to true
    
    for j = 1:numel(orientations)
            
        % PDF
        idx = (trlUncertaintyLevels == i) & (trlStimOris == orientations(j));
        trial_probs(idx) = interp1( ...
            retData.rvOriErrs, ...
            retData.analyticalPDF_stim(j, :), ...
            trlErrors(idx), ...
            'linear');

        % PDF HC
        idxHC = (trlConfReports == 1) & idx;
        trial_probs_HC(idxHC) = interp1( ...
            retData.rvOriErrs, ...
            retData.analyticalPDF_stim_HC(j, :), ...
            trlErrors(idxHC), ...
            'linear');
        trial_probs_Conf(idxHC) = retData.pHC_stim(j);

        % PDF LC
        idxLC = (trlConfReports == 0) & idx;
        trial_probs_LC(idxLC) = interp1( ...
            retData.rvOriErrs, ...
            retData.analyticalPDF_stim_LC(j, :), ...
            trlErrors(idxLC), ...
            'linear');
        trial_probs_Conf(idxLC) = retData.pLC_stim(j); % Mutually exclusive from LC

    end

end

% NLL loss
% P(err|trial) = P(err|trial, confreport)*P(confreport|trial)
ll_HC = log( trial_probs_HC .* trial_probs_Conf + eps); % P(ERR|HC)*P(HC)
ll_LC = log( trial_probs_LC .* trial_probs_Conf + eps); % P(ERR|LC)*P(LC)
% ll    = log( trial_probs + eps ); 

% ll_HC and ll_LC intersection should be zero. Or is it?

% TODO: remove ll - adding this probably leads to overfitting
nll = ( ll_HC + ll_LC); %ll +  % Summation of log => correct thing to do
nll = - sum(nll(:));

end
