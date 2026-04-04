% Loss function for optimization
function nll = computeNLLCov(params, metaData, fitType, optimizationFlag)

if nargin < 4
    optimizationFlag = true;
end

nLevels = metaData.n_levels;

% Params
param_sigma_s        = params(1:nLevels);
param_scale          = params(nLevels + 1);
param_sigma_meta     = params(nLevels + 2);
param_Cc             = params(nLevels + 3);
param_guessrate      = params(nLevels + 4);

if fitType == "full"
    param_ori_scale = params(nLevels + 5);
    param_bias      = params(nLevels + 6);
elseif fitType == "full_jumbo"
    param_a               = params(nLevels + 5:nLevels + 4 + nLevels);
    param_bias            = params(2*nLevels + 5);
    % Hardcoding values for now
    %param_a               = metaData.a; % params(nLevels + 5:nLevels + 4 + nLevels);
    %param_bias            = metaData.biasAmp; %params(2*nLevels + 5);
end

% Metadata
orientations         = metaData.orientations;
trlErrors            = metaData.trlErrors;
trlConfReports       = metaData.trlConfReports;
trlUncertaintyLevels = metaData.trlUncertaintyLevels;
trlStimOris          = metaData.trlStimOris;
sigma_m_shape1       = metaData.sigma_m_shape1;
sigma_m_shape2       = metaData.sigma_m_shape2;
biasShape            = metaData.biasShape;

trial_probs      = zeros(size(trlErrors));
trial_probs_HC   = zeros(size(trlErrors));   % conditional prob of error given conf report - HC
trial_probs_LC   = zeros(size(trlErrors));   % conditional prob of error given conf report - LC
trial_probs_Conf = zeros(size(trlErrors)); % Probability of a trial being high or low confidence

warning("Incorrect nll code. Fix this!")
warning("It might be unfair to compare full vs reduced model.")
warning("To compare these fairly, I need to ensure both models are describing the exact same joint probability space.")

for i=1:nLevels
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = param_guessrate;
    %modelParams.oriErrBinWidth      = 0.1;
    
    if fitType == "full" || fitType == "full_jumbo"

        if fitType == "full"
            modelParams.a              = param_ori_scale.*param_sigma_s(i);
        elseif fitType == "full_jumbo"
            modelParams.a              = param_a(i); 
        end
        
        modelParams.b              = param_sigma_s(i);
        modelParams.biasAmp        = param_bias;
        modelParams.sigma_m_shape1 = sigma_m_shape1;
        modelParams.sigma_m_shape2 = sigma_m_shape2;
        modelParams.biasShape      = biasShape;
        
        retData = getPDFs_cov(orientations, modelParams, true); % Seems like setting this to false makes things slow

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
                'linear'); % this might return zero for out of range values
            trial_probs_Conf(idxLC) = retData.pLC_stim(j); % Mutually exclusive from LC
            % Scale by ori as well

        end

    elseif fitType == "reduced"
        retData = getPDFs_cov_reduced(modelParams, true); % originally set to true
        
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

end

% if fitType == "full" || fitType == "full_jumbo"
%     % P(stim)
%     trial_probs_HC = trial_probs_HC(:)./numel(orientations);
%     trial_probs_LC = trial_probs_LC(:)./numel(orientations);
% end


% NLL loss
% P(err|trial) = P(err|trial, confreport)*P(confreport|trial)
% Very Very Important: USe of eps is problamatic!!!!
% Many small eps can blow up the values. Better otpion might be to clamp
% probabiltites instead
% tmp = trial_probs_HC .* trial_probs_Conf 
% tmp(tmp <= 0) = eps;
% include delta e terms as well ie. in this case 2
% How about conditional NLL?
pHC_ = trial_probs_HC.* trial_probs_Conf;
pHC_ = pHC_(trlConfReports == 1); % this is needed to avoid lot of zeros
pHC_(pHC_ <= 0) = eps;
pLC_ = trial_probs_LC.* trial_probs_Conf;
pLC_ = pLC_(trlConfReports == 0);
pLC_(pLC_ <= 0) = eps;

% I should actually be using much smaller bin size i.e. 0.1
ll_HC = log( pHC_ ); % P(ERR|HC)*P(HC) eps
ll_LC = log( pLC_ ); % P(ERR|LC)*P(LC)
% ll    = log( trial_probs + eps ); 

% ll_HC and ll_LC intersection should be zero. Or is it?

% nll = ( ll_HC + ll_LC); %ll +  % Summation of log => correct thing to do

if ~optimizationFlag
    nll = [-ll_HC' -ll_LC'];
    % nll = - nll(:); % NLL for each trial
else
    % nll = - sum(nll(:)); % Aggregate NLL of all trials
    nll = - sum(ll_HC(:)) - sum(ll_LC(:));
end

end
