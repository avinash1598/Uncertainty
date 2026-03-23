function [retData] = getLLhChoice_MCASANDRE_Contn(modelParams)

    sigma_e                   = modelParams.sigma_e;
    sigma_i                   = modelParams.sigma_i;
    varGain                   = modelParams.varGain;
    internalNoiseSamplesCnt   = modelParams.internalNoiseSamplesCnt;
    Cc                        = modelParams.Cc;
    sigma_m                   = modelParams.sigma_m;
    
    % Get samples of sigma_d
    shapeParam = 1./varGain;
    scaleParam = varGain;
    gammaSamples = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
        shapeParam, scaleParam);
    sigma_theta = sqrt( sigma_e^2 + (sigma_i*gammaSamples).^2 );
    
    % For each value of sigma_theta, find the probability of high
    % confidence and low confidence. These probabilitites are
    % calculated by taking the inverse of the log normal distribution
    % (which basically is a distribution of confidence variable).
    muLogN    = - log((sigma_theta.^2) ./ sqrt(sigma_m.^2 + sigma_theta.^2));
    sigmaLogN = sqrt(log((sigma_m.^2)./(sigma_theta.^2) + 1));

    x1 = repmat(Cc, [1 internalNoiseSamplesCnt]);
    cdf_vals  = logncdf(x1, muLogN, sigmaLogN);
    
    % mean_cdf_val = mean(cdf_vals);
    
    % retData.pHC = 1 - mean_cdf_val;
    % retData.pLC = mean_cdf_val;
    retData.sigma_theta = sigma_theta;
    retData.pHC = 1 - cdf_vals;
    retData.pLC = cdf_vals;
end

