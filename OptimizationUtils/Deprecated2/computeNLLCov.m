% Loss function for optimization
function nll = computeNLLCov(params, metaData, fitType)

nLevels = metaData.n_levels;

% Params
param_sigma_s        = params(1:nLevels);
param_scale          = params(nLevels + 1);
param_sigma_meta     = params(nLevels + 2);
param_Cc             = params(nLevels + 3);
param_guessrate      = params(nLevels + 4);

if fitType == "full"
    %param_ori_scale = params(nLevels + 5);
    %param_bias      = params(nLevels + 6);
    % Hardcoding values for now
    %param_a               = metaData.a; % params(nLevels + 5:nLevels + 4 + nLevels);
    %param_bias            = metaData.biasAmp; %params(2*nLevels + 5);
    param_a               = params(nLevels + 5:nLevels + 4 + nLevels);
    param_bias            = params(2*nLevels + 5);
end

% Metadata
orientations         = metaData.orientations;
trlErrors            = metaData.trlErrors;
trlConfReports       = metaData.trlConfReports;
trlUncertaintyLevels = metaData.trlUncertaintyLevels;
trlStimOris          = metaData.trlStimOris;

if fitType == "full"
    sigma_m_shape1       = metaData.sigma_m_shape1;
    sigma_m_shape2       = metaData.sigma_m_shape2;
    biasShape            = metaData.biasShape;
end    

trial_probs      = zeros(size(trlErrors));
trial_probs_HC   = zeros(size(trlErrors));   % conditional prob of error given conf report - HC
trial_probs_LC   = zeros(size(trlErrors));   % conditional prob of error given conf report - LC
trial_probs_Conf = zeros(size(trlErrors)); % Probability of a trial being high or low confidence

for i=1:nLevels
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = param_guessrate;
    % modelParams.oriErrBinWidth      = 2;
    
    if fitType == "full"
        modelParams.b              = param_sigma_s(i);
        modelParams.a              = param_a(i); %param_ori_scale.*param_sigma_s(i); %param_a(i); 
        modelParams.biasAmp        = param_bias;
        modelParams.sigma_m_shape1 = sigma_m_shape1;
        modelParams.sigma_m_shape2 = sigma_m_shape2;
        modelParams.biasShape      = biasShape;
        
        retData = getEstimationsPDF_cov(orientations, modelParams, true); % Seems like setting this to false makes things slow
        
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

end

% NLL loss

ll_HC = log( trial_probs_HC .* trial_probs_Conf + eps); % P(ERR/HC)*P(HC) or P(ERR/LC)*P(LC)
ll_LC = log( trial_probs_LC .* trial_probs_Conf + eps);
ll    = log( trial_probs + eps ); 

nll = ( ll + ll_HC + ll_LC);
% nll = - sum(nll(:), 'omitnan'); % omitnan or not?? Probably not. It
% messes up with optimization
nll = - sum(nll(:));

end
