function [retData] = getNLL_cov(stimVals, modelParams)

% warning("Deprecated")
if nargin < 4
    optimizationFlag = false;   % default value
end

stimOris                  = stimVals;
b                         = modelParams.b;       % Baseline stimulus dependent sensory noise level
a                         = modelParams.a;       % Amplitude for stimulus dependent sensory noise modulation
biasAmp                   = modelParams.biasAmp;
scale                     = modelParams.scale;
Cc                        = modelParams.Cc;
sigma_meta                = modelParams.sigma_meta;
guessRate                 = modelParams.guessRate;
internalNoiseSamplesCnt   = 500; %1000;
numOris                   = numel(stimOris);

bias = biasAmp*sind(2*stimOris);
% sigma_s_stim = b + a.*(abs(sind(2*stimOris))); %sigma_s_stim = sigma_s_stim';
sigma_s_stim = b + a.*(abs(sind(stimOris - 90))); % as per human subjects

% Internal noise covaries with sensory noise
scaleParam = scale;
shapeParams = sigma_s_stim.^2 ./ scaleParam;
gammaSamples = zeros(numel(shapeParams), internalNoiseSamplesCnt);
% weighted_sigma_m_stim = zeros(numel(shapeParams), internalNoiseSamplesCnt);

for i = 1:numel(shapeParams)
    shapeParam = shapeParams(i);
    gammaSamples(i, :) = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
        shapeParam, scaleParam);

    % weighted_sigma_m_stim(i, :) = gammaSamples(i, :).*gampdf(gammaSamples(i, :), shapeParam, scaleParam);
end

% sigma_m_stim = sqrt( sigma_s_stim'.^2 + (gammaSamples).^2 ); % Measurement noise
sigma_m_stim = sqrt( gammaSamples ); % Measurement noise

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

% Compute PDF for each orientation
analyticalPDF_stim    = zeros(numOris, numel(rvOriErrs));
analyticalPDF_stim_HC = zeros(numOris, numel(rvOriErrs));
analyticalPDF_stim_LC = zeros(numOris, numel(rvOriErrs));

est_sigma_m_stim = zeros(numOris, 1);
est_sigma_m_stim_HC = zeros(numOris, 1);
est_sigma_m_stim_LC = zeros(numOris, 1);

est_mad_m_stim = zeros(numOris, 1);
est_mad_m_stim_HC = zeros(numOris, 1);
est_mad_m_stim_LC = zeros(numOris, 1);

for i = 1:numOris
    [pdf, pdfHC, pdfLC] = getGaussianMixturePDF(rvOriErrs, ...
        sigma_m_stim(i, :), bias(i), pHC_all(i, :), pLC_all(i, :), guessRate);
    
    analyticalPDF_stim(i, :)    = pdf;
    analyticalPDF_stim_HC(i, :) = pdfHC;
    analyticalPDF_stim_LC(i, :) = pdfLC;
    
    est_sigma_m_stim(i) = sigma_;
    est_sigma_m_stim_HC(i) = sigmaHC;
    est_sigma_m_stim_LC(i) = sigmaLC;
    
    est_mad_m_stim(i) = mad_m;
    est_mad_m_stim_HC(i) = mad_m_HC;
    est_mad_m_stim_LC(i) = mad_m_LC;
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

%% Sigma
% Expected sigma - LC and HC by stimulus orientation
% retData.E_sigma_m_stim_HC = sum(sigma_m_stim.*pHC_all, 2)./sum(pHC_all, 2);
% retData.E_sigma_m_stim_LC = sum(sigma_m_stim.*pLC_all, 2)./sum(pLC_all, 2);
retData.E_sigma_m_stim_HC = est_sigma_m_stim_HC;
retData.E_sigma_m_stim_LC = est_sigma_m_stim_LC;

% Expected sigma - LC and HC - combined across all orientations
% retData.E_sigma_m_HC = sum(sigma_m_stim.*pHC_all, 'all')/sum(pHC_all, 'all');
% retData.E_sigma_m_LC = sum(sigma_m_stim.*pLC_all, 'all')/sum(pLC_all, 'all');
retData.E_sigma_m_HC = sqrt( sum( est_sigma_m_stim_HC.^2.*retData.pHC_stim )./sum(retData.pHC_stim) );
retData.E_sigma_m_LC = sqrt( sum( est_sigma_m_stim_LC.^2.*retData.pLC_stim )./sum(retData.pLC_stim) );

% % Expected sigma for each stimulus orientation - aggreagted by HC and LC:
% Analytical solution
% retData.E_sigma_m_stim = sigma_s_stim;
retData.E_sigma_m_stim = est_sigma_m_stim;
% retData.E_sigma_m_stim = sqrt( scaleParam^2 .* sigma_s_stim.^2 .* (sigma_s_stim.^2 / scaleParam + 1) );

% Expected sigma for all perceptual reports
% retData.E_sigma_m = sqrt( mean( sigma_s_stim.^2 ) );
retData.E_sigma_m = sqrt( mean( est_sigma_m_stim.^2 ) );% + std(bias).^2 ;  
% retData.E_sigma_m = sqrt( mean( scaleParam^2 .* sigma_s_stim.^2 .* (sigma_s_stim.^2 / scaleParam + 1) ) );

%% MAD
retData.mad_m_stim_HC = est_mad_m_stim_HC;
retData.mad_m_stim_LC = est_mad_m_stim_LC;
retData.mad_m_stim    = est_mad_m_stim;

retData.mad_m_HC     = ( sum( est_mad_m_stim_HC.*retData.pHC_stim )./sum(retData.pHC_stim) );
retData.mad_m_LC     = ( sum( est_mad_m_stim_LC.*retData.pLC_stim )./sum(retData.pLC_stim) );
retData.mad_m        = mean( est_mad_m_stim); % + mad(bias)
retData.mad_m_by_ori = est_mad_m_stim; % change 1

%%
% Stimulus dependent bias
% retData.bias = biasAmp*sind(2*stimOris);
retData.bias = bias; %biasAmp*sind(2*stimOris);

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
function [pdf, pdfHC, pdfLC] = getGaussianMixturePDF( ...
    x, sigma_m_stim, bias, pHC, pLC, guessRate)
% x = array of perceptual errors

pdf_random_guesses = ones(size(x)) / (numel(x) + 1);
pdf_random_guesses = pdf_random_guesses./trapz(x, pdf_random_guesses); 
pdf_random_guesses = pdf_random_guesses';

% I guess using trapz is fine. Error reports lives in that range
% Note: wrapped distribution is conceptually correct - but if the error is
% not too large then truncated distribution is fine 
% (truncated might not be best for the highest uncertainty level)
Z = x - bias;
p_X = exp(- (Z').^2 ./ (2*sigma_m_stim.^2 ) ) ./ sqrt(2*pi*sigma_m_stim.^2);
pdf = sum(p_X, 2);
pdf = pdf./trapz(x, pdf);
pdf = (1 - guessRate)*pdf + guessRate*pdf_random_guesses;
pdf = pdf./trapz(x, pdf);

% HC pdf
p_X_HC = p_X.*pHC;
pdfHC = sum(p_X_HC, 2);
pdfHC = pdfHC./trapz(x, pdfHC);
% pdfHC = (1 - guessRate)*pdfHC + guessRate*pdf_random_guesses; % (No guess rate for high confidence)
% pdfHC = pdfHC./trapz(x, pdfHC);

% LC pdf
p_X_LC = p_X.*pLC;
pdfLC = sum(p_X_LC, 2);
pdfLC = pdfLC./trapz(x, pdfLC);
pdfLC = (1 - guessRate)*pdfLC + guessRate*pdf_random_guesses; % Guess rate for low confidence
pdfLC = pdfLC./trapz(x, pdfLC);

end