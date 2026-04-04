% Loss function for optimization
function nll = computeNLLCov_binned(params, metaData, fitType, optimizationFlag)

if nargin < 4
    optimizationFlag = true;
end

nLevels = metaData.n_levels;

% Params
param_sigma_s        = params(1:nLevels);
param_scale          = params(nLevels + 1);
param_sigma_meta     = params(nLevels + 2);
param_Cc             = params(nLevels + 3);
param_exp            = params(nLevels + 4);
param_guessrate      = params(nLevels + 5);

if fitType == "full"
    param_ori_scale = params(nLevels + 6);
    param_bias      = params(nLevels + 7);
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

totalNLL = 0;
trlNLLs = zeros(size(trlErrors));

% trial_probs      = zeros(size(trlErrors));
% trial_probs_HC   = zeros(size(trlErrors));   % conditional prob of error given conf report - HC
% trial_probs_LC   = zeros(size(trlErrors));   % conditional prob of error given conf report - LC
% trial_probs_Conf = zeros(size(trlErrors)); % Probability of a trial being high or low confidence

warning("Incorrect nll code. Fix this!")
warning("It might be unfair to compare full vs reduced model.")
warning("To compare these fairly, I need to ensure both models are describing the exact same joint probability space.")

for i=1:nLevels
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = param_guessrate;
    modelParams.exp                 = param_exp;
    %modelParams.oriErrBinWidth      = 0.1;
    
    if fitType == "full"

        modelParams.a              = param_ori_scale.*param_sigma_s(i);
        modelParams.b              = param_sigma_s(i);
        modelParams.biasAmp        = param_bias;
        modelParams.sigma_m_shape1 = sigma_m_shape1;
        modelParams.sigma_m_shape2 = sigma_m_shape2;
        modelParams.biasShape      = biasShape;
        
        retData = getPDFs_cov(orientations, modelParams, true); % Seems like setting this to false makes things slow

        for j = 1:numel(orientations)
            
            idx = (trlUncertaintyLevels == i) & (trlStimOris == orientations(j));
             
            % PDF HC
            idxHC = (trlConfReports == 1) & idx;
            trlProbsHC = interp1( ...
                retData.rvOriErrs, ...
                retData.analyticalPDF_stim_HC(j, :), ...
                trlErrors(idxHC), ...
                'linear');
            probHC = retData.pHC_stim(j);
            tmp = trlProbsHC*probHC; tmp(tmp <= 0) = eps;
            nllsHC = -log( tmp );

            % PDF LC
            idxLC = (trlConfReports == 0) & idx;
            trlProbsLC = interp1( ...
                retData.rvOriErrs, ...
                retData.analyticalPDF_stim_LC(j, :), ...
                trlErrors(idxLC), ...
                'linear'); % this might return zero for out of range values
            probLC = retData.pLC_stim(j); % Mutually exclusive from LC
            tmp = trlProbsLC*probLC; tmp(tmp <= 0) = eps;
            nllsLC = -log( tmp );
            
            totalNLL = totalNLL + sum(nllsHC) + sum(nllsLC);
            trlNLLs(idxHC) = nllsHC;
            trlNLLs(idxLC) = nllsLC;
            
            % Scale by ori as well - Not needed

        end
    end

end

if ~optimizationFlag
    nll = trlNLLs;
else
    nll = totalNLL;
end

end
