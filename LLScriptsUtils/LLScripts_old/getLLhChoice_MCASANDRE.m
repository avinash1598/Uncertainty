function [retData] = getLLhChoice_MCASANDRE(stimVals, modelParams)

    sigma_e                   = modelParams.sigma_e;
    sigma_i                   = modelParams.sigma_i;
    varGain                   = modelParams.varGain;
    internalNoiseSamplesCnt   = modelParams.internalNoiseSamplesCnt;
    Cd                        = modelParams.Cd;
    Cc                        = modelParams.Cc;
    sigma_m                   = modelParams.sigma_m;
    sampleCnt                 = modelParams.sampleCnt;
    
    choicePDFs = zeros(4, numel(stimVals));
    
    for sidx = 1:numel(stimVals)
        stimVal = stimVals(sidx);
        
        % Get samples of sigma_d
        shapeParam = 1./varGain;
        scaleParam = varGain;
        gammaSamples = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
            shapeParam, scaleParam);
        sigma_d = sqrt( sigma_e^2 + (sigma_i*gammaSamples).^2 );
        
        cdf_vals = zeros(internalNoiseSamplesCnt, sampleCnt, 3);
        
        for sdIdx=1:numel(sigma_d)
        
            % Uniformly sample lognormal cdf to obtain estimated of sigma_d evenly
            % tiling the distribution. TODO: try weighted averaging but probabilities
            % around each sample is same, so this step can just be omitted
            muLogN    = log((sigma_d(sdIdx).^2) ./ sqrt(sigma_m.^2 + sigma_d(sdIdx).^2));
            sigmaLogN = sqrt(log((sigma_m.^2)./(sigma_d(sdIdx).^2) + 1));
            sigma_d_hat_samples  = logninv(linspace(1/sampleCnt, 1 - 1/sampleCnt, sampleCnt), muLogN, sigmaLogN);
            
            % Mean and var of confidence variable
            mu_c = (stimVal - Cd) ./ sigma_d_hat_samples;
            sigma_c = sigma_d(sdIdx) ./ sigma_d_hat_samples;
            
            % These three CDF values will be used to compute the probabilitites
            x = [-Cc, 0, Cc];
            x1 = repmat(x, [sampleCnt 1]);
            x2 = repmat(mu_c', [1 numel(x)]);
            x3 = repmat(sigma_c', [1 numel(x)]);
            cdf_vals(sdIdx, :, :) = normcdf(x1, x2, x3);
        end
        
        % TODO: weight cdf_vals by corresponding probabilitites. This probably
        % can be omitted since probabilitites around each sample is the same.
        mean_cdf_vals = mean( squeeze(mean(cdf_vals, 1)), 1 );
        mean_cdf_vals = mean(mean_cdf_vals, 1);
        
        p_c1_d0 = mean_cdf_vals(1);
        p_c0_d0 = mean_cdf_vals(2) - mean_cdf_vals(1);
        p_c1_d1 = 1 - mean_cdf_vals(3);
        p_c0_d1 = mean_cdf_vals(3) - mean_cdf_vals(2);
    
        choicePDFs(1, sidx) = p_c1_d0;
        choicePDFs(2, sidx) = p_c0_d0;
        choicePDFs(3, sidx) = p_c1_d1;
        choicePDFs(4, sidx) = p_c0_d1;
    
    end
    
    retData = choicePDFs;
end

