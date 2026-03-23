function [retData] = getLLh_EstimationTask(rvOriErrs, stimVals, modelParams)

% rvOriErrs - random variable for orientations errors over which
% PDF/likelihood is computed.

stimOris                  = stimVals;
b                         = modelParams.b;       % Baseline stimulus dependent sensory noise level
a                         = modelParams.a;       % Amplitude for stimulus dependent sensory noise modulation
biasAmp                   = modelParams.biasAmp;
shape                     = modelParams.shape;
scale                     = modelParams.scale;
Cc                        = modelParams.Cc;
sigma_meta                = modelParams.sigma_meta;
internalNoiseSamplesCnt   = 1000;
numOris                   = numel(stimOris);  

bias = biasAmp*sind(2*stimOris);

shapeParam = shape;
scaleParam = scale;
gammaSamples = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
    shapeParam, scaleParam);

sigma_s_stim = b + a.*(abs(sind(2*stimOris))); %sigma_s_stim = sigma_s_stim';
sigma_m_stim = sqrt( sigma_s_stim'.^2 + (gammaSamples).^2 ); % Measurement noise

% For each value of sigma_m_stim, find the probability of high
% confidence and low confidence. These probabilitites are
% calculated by taking the inverse of the log normal distribution
% (which basically is a distribution of confidence variable).
muLogN    = - log((sigma_m_stim.^2) ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2));
sigmaLogN = sqrt(log((sigma_meta.^2)./(sigma_m_stim.^2) + 1));

x1 = repmat(Cc, [numOris internalNoiseSamplesCnt]);
cdf_vals  = logncdf(x1, muLogN, sigmaLogN);

pHC_all = 1 - cdf_vals;
pLC_all = cdf_vals;

% Calculate probability of HC and LC split by orientation. 
% While average should I consider PDF
% values, since the chances of a particular gammasamples are not equal? No!
% LLN is the key.
mean_cdf_val = mean(cdf_vals(:));
mean_cdf_val_stim = mean(cdf_vals, 2);

% Weight prob HC and LC by prob of corresponding sample
% Note: result of this will literally be same as pLC_byOri when no 
% of gammasamples are large. 
pdf = gampdf(gammaSamples, shapeParam, scaleParam);
y = cdf_vals.*pdf;
I1 = trapz(gammaSamples, y, 2); % Integrate y with non-uniform spacing gammaSamples
% I2 = trapz(stimOris, I1*(1/numOris));

% Is I2 correct way to find the expectation of PLC? I2 might be equal to
% mean_cdf_val only when no of thetas is large. Run simulations to see what
% is analytical solution. I would assume mean_cdf_val to be the correct solution.

diff = abs(mean_cdf_val_stim - I1);
assert(all(diff < 0.1), 'Sanity check failed: not all differences < 0.1'); % Make sure the two values are same (just a sanity check for WLLN)

% Compute PDF for each orientation
analyticalPDF_stim    = zeros(numOris, numel(rvOriErrs));
analyticalPDF_stim_HC = zeros(numOris, numel(rvOriErrs));
analyticalPDF_stim_LC = zeros(numOris, numel(rvOriErrs));

for i = 1:numOris
    [pdf, pdfHC, pdfLC] = getGaussianMixturePDF(rvOriErrs, ...
        sigma_m_stim(i, :), bias(i), pHC_all(i, :), pLC_all(i, :));

    analyticalPDF_stim(i, :)    = pdf;
    analyticalPDF_stim_HC(i, :) = pdfHC;
    analyticalPDF_stim_LC(i, :) = pdfLC;
end

% Compute PDF for all orientation
analyticalPDF = mean(analyticalPDF_stim, 1);
analyticalPDF = analyticalPDF / trapz(rvOriErrs, analyticalPDF);

% PDF for HC and LC
analyticalPDF_HC = mean(analyticalPDF_stim_HC, 1);
analyticalPDF_HC = analyticalPDF_HC / trapz(rvOriErrs, analyticalPDF_HC);

analyticalPDF_LC = mean(analyticalPDF_stim_LC, 1);
analyticalPDF_LC = analyticalPDF_LC / trapz(rvOriErrs, analyticalPDF_LC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Return data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Probability of HC and LC (Split by stimulus orientation)
retData.pHC_stim   = 1 - mean_cdf_val_stim;
retData.pLC_stim   = mean_cdf_val_stim;

% Probability of HC and LC (combined across all orientations)
retData.pHC  = 1 - mean_cdf_val;
retData.pLC  = mean_cdf_val;

% Expected sigma - LC and HC by stimulus orientation
retData.E_sigma_m_stim_HC = sum(sigma_m_stim.*pHC_all, 2)./sum(pHC_all, 2);
retData.E_sigma_m_stim_LC = sum(sigma_m_stim.*pLC_all, 2)./sum(pLC_all, 2);

% Expected sigma - LC and HC - combined across all orientations
retData.E_sigma_m_HC = sum(sigma_m_stim.*pHC_all, 'all')/sum(pHC_all, 'all');
retData.E_sigma_m_LC = sum(sigma_m_stim.*pLC_all, 'all')/sum(pLC_all, 'all');

% Expected sigma for each stimulus orientation - aggreagted by HC and LC:
% Analytical solution
retData.E_sigma_m_stim = sqrt( sigma_s_stim.^2 + scaleParam.^2 * (shapeParam * (shapeParam + 1)) );

% Expected sigma for all perceptual reports
% retData.E_sigma_m = sqrt( mean( sigma_s_stim.^2 ) + sigma_si.^2 * scaleParam.^2 * (shapeParam * (shapeParam + 1)) );
retData.E_sigma_m = sqrt( mean( sigma_s_stim.^2 + scaleParam.^2 * (shapeParam * (shapeParam + 1)) ) );

% Stimulus dependent bias
% retData.bias = biasAmp*sind(4*stimOris);
retData.bias = bias;

% PDF for each orientation
retData.rvOriErrs = rvOriErrs;

retData.analyticalPDF_stim = analyticalPDF_stim;
retData.analyticalPDF_stim_LC = analyticalPDF_stim_LC; % LC
retData.analyticalPDF_stim_HC = analyticalPDF_stim_HC; % HC

% PDF combined across all orientations
retData.analyticalPDF = analyticalPDF;
retData.analyticalPDF_LC = analyticalPDF_LC; % LC
retData.analyticalPDF_HC = analyticalPDF_HC; % HC

end

% PDF for individual orientations
function [pdf, pdfHC, pdfLC] = getGaussianMixturePDF(x, sigma_m_stim, bias, pHC, pLC)
% x = array of perceptual errors

x = x - bias;
p_X = exp(- (x').^2 ./ (2*sigma_m_stim.^2 ) ) ./ sqrt(2*pi*sigma_m_stim.^2);
pdf = sum(p_X, 2);
pdf = pdf./trapz(x, pdf);

% HC pdf
p_X_HC = p_X.*pHC;
pdfHC = sum(p_X_HC, 2);
pdfHC = pdfHC./trapz(x, pdfHC);

% LC pdf
p_X_LC = p_X.*pLC;
pdfLC = sum(p_X_LC, 2);
pdfLC = pdfLC./trapz(x, pdfLC);

end

