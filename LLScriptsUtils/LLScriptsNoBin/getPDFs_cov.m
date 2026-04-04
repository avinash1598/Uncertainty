function [retData] = getPDFs_cov(stimVals, modelParams, optimizationFlag)

% warning("Deprecated. Use LLScriptsNoBin instead.")

% Commenting it out to save execution time
if nargin < 3
    optimizationFlag = false;   % default value
end

% Code can be made fast by omitting if else conditions
% No of samples to draw for numerical approximation
if ~isfield(modelParams, 'internalNoiseSamplesCnt')
    modelParams.internalNoiseSamplesCnt = 100; % Default value obtained from params check analysis
end

% Precision of orientation error bins
if ~isfield(modelParams, 'oriErrBinWidth')
    modelParams.oriErrBinWidth = 2; %2; % Default value obtained from params check analysis
end

% Controls wrapped gaussian distribution
if ~isfield(modelParams, 'K')
    modelParams.K = 1; % Default value obtained from params check analysis
end

K                         = modelParams.K;
rvOriErrs                 = -90:modelParams.oriErrBinWidth:90; % imp to generate bw -90 and 90 to avoid NLL issues
stimOris                  = stimVals;
b                         = modelParams.b;       % Baseline stimulus dependent sensory noise level
a                         = modelParams.a;       % Amplitude for stimulus dependent sensory noise modulation
biasAmp                   = modelParams.biasAmp;
scale                     = modelParams.scale;
Cc                        = modelParams.Cc;
sigma_meta                = modelParams.sigma_meta;
guessRate                 = modelParams.guessRate;
internalNoiseSamplesCnt   = modelParams.internalNoiseSamplesCnt; % 100 values is good enough (irrespctive of std of data)
numOris                   = numel(stimOris);
sigma_m_shape1            = modelParams.sigma_m_shape1;
sigma_m_shape2            = modelParams.sigma_m_shape2;
biasShape                 = modelParams.biasShape;
nonlinearityExp           = modelParams.exp;

bias         = biasAmp*sind(biasShape*stimOris);
sigma_s_stim = b + a.*(abs(sind(sigma_m_shape1*stimOris - sigma_m_shape2)));
% bias = biasAmp*sind(2*stimOris);
% sigma_s_stim = b + a.*(abs(sind(stimOris - 90))); % as per human subjects
% % sigma_s_stim = b + a.*(abs(sind(2*stimOris))); %sigma_s_stim = sigma_s_stim';

% % Internal noise covaries with sensory noise
% scaleParam = scale;
% shapeParams = sigma_s_stim.^2 ./ scaleParam;
% gammaSamples = zeros(numel(shapeParams), internalNoiseSamplesCnt);

% for i = 1:numel(shapeParams)
%     shapeParam = shapeParams(i);
%     gammaSamples(i, :) = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
%         shapeParam, scaleParam);
% end
% 
% sigma_m_stim = sqrt( gammaSamples ); % Measurement noise

% Sample from log normal distribution instead
lognSamples = zeros(numOris, internalNoiseSamplesCnt);
u = linspace(1/internalNoiseSamplesCnt, 1-1/internalNoiseSamplesCnt, internalNoiseSamplesCnt);

for i = 1:numOris
    % Find parameters of log normal distribution (i.e. mean and variance of normal distribution I guess)
    % This parameter is what matlab accepts
    % mu and sigma here are the parameters and not necesarily the actual mean
    % and variance
    m = sigma_s_stim(i).^2; v = scale*(m^nonlinearityExp); % desired mean & variance of lognormal distribution
    sigma2 = log(1 + v/m^2);  % variance
    mu = log(m) - 0.5*sigma2; %mean
    
    lognSamples(i, :) = logninv(u, mu, sqrt(sigma2));
end

sigma_m_stim = sqrt( lognSamples ); % Measurement noise

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
    [pdf, pdfHC, pdfLC, sigma_, sigmaHC, sigmaLC, mad_m, mad_m_HC, mad_m_LC] = getGaussianMixturePDF(rvOriErrs, ...
        sigma_m_stim(i, :), bias(i), pHC_all(i, :), pLC_all(i, :), guessRate, ...
        K, ...
        optimizationFlag);

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

% PDF for HC and LC
analyticalPDF_HC = mean(analyticalPDF_stim_HC, 1);
analyticalPDF_LC = mean(analyticalPDF_stim_LC, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Return data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Probability of HC and LC (Split by stimulus orientation)
retData.pHC_stim   = 1 - mean_cdf_val_stim;
retData.pLC_stim   = mean_cdf_val_stim;

% Probability of HC and LC (combined across all orientations)
retData.pHC  = 1 - mean_cdf_val;
retData.pLC  = mean_cdf_val;

if ~optimizationFlag
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
end

%%
% Stimulus dependent bias
% retData.bias = biasAmp*sind(2*stimOris);

% PDF for each orientation
retData.rvOriErrs = rvOriErrs;

retData.analyticalPDF_stim = analyticalPDF_stim;
retData.analyticalPDF_stim_LC = analyticalPDF_stim_LC; % LC
retData.analyticalPDF_stim_HC = analyticalPDF_stim_HC; % HC

% PDF combined across all orientations
retData.analyticalPDF = analyticalPDF;
retData.analyticalPDF_LC = analyticalPDF_LC; % LC
retData.analyticalPDF_HC = analyticalPDF_HC; % HC

retData.analytical_sigma_s_stim = sigma_s_stim;
retData.bias = bias; 
retData.analytical_bias = bias;
retData.analytical_sigma_s_reduced = sqrt( mean( sigma_s_stim.^2, 2 ) + std(bias).^2 );

% this is sending a lot of information
end

% PDF for individual orientations
function [pdf, pdfHC, pdfLC, sigma, sigmaHC, sigmaLC, mad_m, mad_m_HC, mad_m_LC] = getGaussianMixturePDF( ...
    x, sigma_m_stim, bias, pHC, pLC, guessRate, K, optimizationFlag)
% x = array of perceptual errors

% If probability lands outside [-90,90] i.e. period 180, move it back inside by wrapping.
period = 180;
% K = 1; % K=1 good enough for max obsevred signam of 45 in the actual data. If needed maybe 2 can be tried as well
pdf   = zeros(length(x),1);
pdfHC = zeros(length(x),1);
pdfLC = zeros(length(x),1);

% Normalize such that sum is 1
% This is done to avoid using trapz
% keyboard

% Two properties needs to be satisfied
% 1. pHC + pLC = const for all ideally 1
% 2. sum(pHC) and sum(PLC) = 1
pHC_norm = pHC./sum(pHC); 
pLC_norm = pLC./sum(pLC);
% norm_const = (sum(pHC) + sum(pLC));
% pHC_norm = pHC./norm_const; 
% pLC_norm = pLC./norm_const;

inv2sig2  = 1 ./ (2*sigma_m_stim.^2);
normconst = 1 ./ (sqrt(2*pi) .* sigma_m_stim);

for k = -K:K
    Z = x - bias + k*period;
    p_X = exp(-(Z'.^2).*inv2sig2) .* normconst;
    % pdf = pdf + sum(p_X,2);
    pdf = pdf + mean(p_X,2);
    
    % HC
    % This seems some sort of weighted sum maybe!
    p_X_HC = p_X.*pHC_norm;
    pdfHC  = pdfHC + sum(p_X_HC, 2);
    
    % LC
    p_X_LC = p_X.*pLC_norm;
    pdfLC  = pdfLC + sum(p_X_LC, 2);
    
end

% PDF taking guess rate into account
% Effective guess rate taking probability of low confidence into account
pdfRG = 1/180;
guessRateEff = guessRate*mean(pLC);
pdf = (1 - guessRateEff)*pdf + guessRateEff*pdfRG;

% PDF LC
pdfRG = 1/180;
pdfLC = (1 - guessRate)*pdfLC + guessRate*pdfRG;

% pdfLC = pdfLC./trapz(x,pdfLC);
% pdfHC = pdfHC./trapz(x,pdfHC);

% keyboard
% disp(trapz(x,pdf))
% disp(trapz(x,pdfHC))
% disp(trapz(x,pdfLC))

areaHC = trapz(x,pdfHC);
areaLC = trapz(x,pdfLC);
tol = 1e-1;

assert( abs(areaHC - 1) < tol );
assert( abs(areaLC - 1) < tol );

% assert to make sure it is very close to 1;

%% MAD and STD
% Calculate metric only after optimization is complete
% Compute wrt to zero (assume zero to be mean and the median)
if ~optimizationFlag
    
    %% Calculate std dev
    dx = x(2) - x(1); % Assuming uniform
    % warning("Make sure error bins are uniformly spaced")
    % Compute wrt to zero to account for bias
    mu = 0; %sum( x.*pdf'*dx );
    sigma = sqrt( sum( ((x - mu).^2).*pdf'*dx ) );
    
    % For HC
    mu = 0;%sum( x.*pdfHC'*dx );
    sigmaHC = sqrt( sum( ((x - mu).^2).*pdfHC'*dx ) );
    
    % For LC
    mu = 0; %sum( x.*pdfLC'*dx );
    sigmaLC = sqrt( sum( ((x - mu).^2).*pdfLC'*dx ) );

    %% MAD from PDF
    dx = x(2)-x(1);
    F = cumsum(pdf) * dx;
    if ~isnan(F(end))
        F = F / F(end);   % normalize
        % [Funiq, idx] = unique(F);
        % xuniq = x(idx);
        median_val = 0; %interp1(Funiq, xuniq, 0.5);
        mad_fun = @(d) (interp1(x, F, median_val + d) - interp1(x, F, median_val - d)) - 0.5;
        
        d0 = (x(end) - x(1)) / 4; % initial guess for d
        MAD = fzero(mad_fun, d0);
        mad_m = MAD;
    else
        mad_m = nan;
    end
    
    % HC
    dx = x(2)-x(1);
    F = cumsum(pdfHC) * dx;
    if ~isnan(F(end))
        F = F / F(end);   % normalize
        % [Funiq, idx] = unique(F);
        % xuniq = x(idx);
        median_val = 0; %interp1(Funiq, xuniq, 0.5);
        mad_fun = @(d) (interp1(x, F, median_val + d) - interp1(x, F, median_val - d)) - 0.5;
        
        d0 = (x(end) - x(1)) / 4; % initial guess for d
        MAD = fzero(mad_fun, d0);
        mad_m_HC = MAD;
    
    else
        mad_m_HC = nan; % dummy value to avoid nan error
    end
    
    % LC
    dx = x(2)-x(1);
    F = cumsum(pdfLC) * dx;
    if ~isnan(F(end))
        F = F / F(end);   % normalize
        % [Funiq, idx] = unique(F);
        % xuniq = x(idx);
        median_val = 0; %interp1(Funiq, xuniq, 0.5); 
        mad_fun = @(d) (interp1(x, F, median_val + d) - interp1(x, F, median_val - d)) - 0.5;
        
        d0 = (x(end) - x(1)) / 4; % initial guess for d
        MAD = fzero(mad_fun, d0);
        mad_m_LC = MAD;
    else
        mad_m_LC = nan;
    end

else
    mad_m = 0;
    mad_m_HC = 0;
    mad_m_LC = 0;
    sigma = 0;
    sigmaHC = 0;
    sigmaLC = 0;
end

end
