% Loss function for optimization
function nll = computeNLL(params, metaData, fitType)

nLevels = metaData.n_levels;

% Params
param_sigma_s         = params(1:nLevels);
param_shape           = params(nLevels + 1);
param_scale           = params(nLevels + 2);
param_sigma_meta      = params(nLevels + 3);
param_Cc              = params(nLevels + 4);
param_guessrate       = params(nLevels + 5);

if fitType == "full"
    param_sigma_ori_scale = params(nLevels + 6);
    param_bias            = params(nLevels + 7);
end

% Metadata
errBins        = metaData.errBins;
binned_err_HC  = metaData.binned_err_HC;
binned_err_LC  = metaData.binned_err_LC;
targetMADs     = metaData.targetMADs;
targetMADs_HC  = metaData.targetMADs_HC;
targetMADS_LC  = metaData.targetMADS_LC;
hyperParamC1   = metaData.hyperParamC1;

currPdfFit_HC = zeros(nLevels, numel(errBins));
currPdfFit_LC = zeros(nLevels, numel(errBins));
curr_pHC       = zeros(nLevels, 1);
curr_pLC       = zeros(nLevels, 1);
curr_mad_m      = zeros(nLevels, 1);
curr_mad_m_HC   = zeros(nLevels, 1);
curr_mad_m_LC   = zeros(nLevels, 1);

for i=1:nLevels
    
    modelParams.sigma_s             = param_sigma_s(i);
    modelParams.shape               = param_shape;   
    modelParams.scale               = param_scale;
    modelParams.Cc                  = param_Cc;
    modelParams.sigma_meta          = param_sigma_meta;
    modelParams.guessRate           = param_guessrate;
    
    if fitType == "full"
        modelParams.b             = param_sigma_s(i);
        modelParams.a             = param_sigma_ori_scale*param_sigma_s(i);
        modelParams.biasAmp       = param_bias;

        retData = getEstimatesPDFs(1:10:180, errBins, modelParams, false);
    else
        % retData = getEstimatesPDFs_reduced_model(errBins, modelParams, false); % originally set to true
    end
    
    % Data for NLL
    currPdfFit_HC(i, :) = retData.analyticalPDF_HC;
    currPdfFit_LC(i, :) = retData.analyticalPDF_LC;
    curr_pHC(i)         = retData.pHC;
    curr_pLC(i)         = retData.pLC;
    curr_mad_m(i)       = retData.mad_m;
    curr_mad_m_HC(i)    = retData.mad_m_HC;
    curr_mad_m_LC(i)    = retData.mad_m_LC;
end

constraint = sum( ( curr_mad_m - targetMADs ).^2 ) + ...
    sum( ( curr_mad_m_HC - targetMADs_HC ).^2 ) + ...
    sum( ( curr_mad_m_LC - targetMADS_LC ).^2 );

% NLL loss
ll_HC = binned_err_HC .* log( currPdfFit_HC.*curr_pHC + eps );
ll_LC = binned_err_LC .* log( currPdfFit_LC.*curr_pLC + eps ); 

nll = - ( sum(ll_HC(:)) + sum(ll_LC(:)) ) + hyperParamC1*constraint;

% fprintf('\n%.2f, %.2f \n',  - ( sum(ll_HC(:)) + sum(ll_LC(:)) ), hyperParamC1*constraint);

end