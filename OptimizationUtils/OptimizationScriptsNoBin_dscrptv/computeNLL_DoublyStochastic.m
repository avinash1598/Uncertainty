% Loss function for optimization
function nll = computeNLL_DoublyStochastic(params, metaData, fitType, optimizationFlag)

if nargin < 4
    optimizationFlag = true;
end

nLevels = metaData.n_levels;

% Params
param_sigma_s        = params(1);
param_scale          = params(2);
param_sigma_meta     = params(3);
param_Cc             = params(4);
param_guessrate      = params(5);

if fitType == "full"
    param_ori_scale = params(6);
    param_bias      = params(7);
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

totalNLL      = 0;
trialLLs      = zeros(size(trlErrors));

% warning("Incorrect nll code. Fix this!")
% warning("It might be unfair to compare full vs reduced model.")
% warning("To compare these fairly, I need to ensure both models are describing the exact same joint probability space.")

for i=1:nLevels
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = param_guessrate;
    
    if fitType == "full"

        modelParams.a              = param_ori_scale.*param_sigma_s;
        modelParams.b              = param_sigma_s;
        modelParams.biasAmp        = param_bias;
        modelParams.sigma_m_shape1 = sigma_m_shape1;
        modelParams.sigma_m_shape2 = sigma_m_shape2;
        modelParams.biasShape      = biasShape;
        
        trlErrs_L      = trlErrors(trlUncertaintyLevels == i);
        trlConf_L      = trlConfReports(trlUncertaintyLevels == i);
        trlStimOri_L   = trlStimOris(trlUncertaintyLevels == i);
        
        retData = getLLs_DoublyStochastic(trlErrs_L, ...
            trlConf_L, ...
            trlStimOri_L, ...
            orientations, ...
            modelParams); % Seems like setting this to false makes things slow
        
        totalNLL = totalNLL + retData.totalNLL;
        trialLLs(trlUncertaintyLevels == i) = retData.trlNLLs;
    end

end


if ~optimizationFlag
    nll = trialLLs;
else
    nll = totalNLL;
end

end
