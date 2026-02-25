function [retData] = getEstimationsPDF_cov(stimVals, rvOriErrs, modelParams, optimizationFlag)

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
internalNoiseSamplesCnt   = 1000;
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

% N = internalNoiseSamplesCnt;
% M = numel(shapeParams);
% 
% p = linspace(1/N, 1-1/N, N);
% P = repmat(p, M, 1);
% A = repmat(shapeParams(:), 1, N);
% 
% gammaSamples = gaminv(P, A, scaleParam);

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

% [pdf_s, pdfHC_s, pdfLC_s, sigma_s, sigmaHC_s, sigmaLC_s, mad_m_s, mad_m_HC_s, mad_m_LC_s] = getGaussianMixturePDF( ...
%     rvOriErrs, sigma_m_stim, bias, pHC_all, pLC_all, guessRate, optimizationFlag);
% 
% analyticalPDF_stim(:)    = pdf_s;
% analyticalPDF_stim_HC(:) = pdfHC_s;
% analyticalPDF_stim_LC(:) = pdfLC_s;
% 
% est_sigma_m_stim(:)    = sigma_s;
% est_sigma_m_stim_HC(:) = sigmaHC_s;
% est_sigma_m_stim_LC(:) = sigmaLC_s;
% 
% est_mad_m_stim(:)    = mad_m_s;
% est_mad_m_stim_HC(:) = mad_m_HC_s;
% est_mad_m_stim_LC(:) = mad_m_LC_s;

for i = 1:numOris
    [pdf, pdfHC, pdfLC, sigma_, sigmaHC, sigmaLC, mad_m, mad_m_HC, mad_m_LC] = getGaussianMixturePDF(rvOriErrs, ...
        sigma_m_stim(i, :), bias(i), pHC_all(i, :), pLC_all(i, :), guessRate, optimizationFlag);
    
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
% retData.bias = biasAmp*sind(4*stimOris);
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


% % PDF for all orientations
% function [pdf_s, pdfHC_s, pdfLC_s, sigma_s, sigmaHC_s, sigmaLC_s, mad_m_s, mad_m_HC_s, mad_m_LC_s] = getGaussianMixturePDF( ...
%     x, sigma_m_stim, bias, pHC, pLC, guessRate, optimizationFlag)
% % x = array of perceptual errors
% 
% pdf_random_guesses = ones(size(x)) / (numel(x) + 1);
% pdf_random_guesses = pdf_random_guesses./trapz(x, pdf_random_guesses); 
% % pdf_random_guesses = pdf_random_guesses';
% 
% % tic
% Z = x - bias';
% A = reshape(Z, size(sigma_m_stim, 1), numel(x), 1);
% B = reshape(sigma_m_stim, size(sigma_m_stim, 1), 1, size(sigma_m_stim, 2));
% p_X = exp(- (A).^2 ./ (2*B.^2 ) ) ./ sqrt(2*pi*B.^2);
% % elapsed_time = toc;
% % disp(['Execution time::: ', num2str(elapsed_time), ' seconds']);
% 
% pdf = sum(p_X, 3);
% area = trapz(x, pdf, 2);
% pdf = pdf./area;
% pdf_s = (1 - guessRate)*pdf + guessRate*pdf_random_guesses; 
% % pdf = pdf./trapz(x, pdf);
% 
% % HC pdf
% A = reshape(pHC, size(pHC, 1), 1, size(pHC, 2));
% p_X_HC = p_X.*A;
% pdfHC = sum(p_X_HC, 3);
% pdfHC = pdfHC./trapz(x, pdfHC, 2);
% pdfHC_s = (1 - guessRate)*pdfHC + guessRate*pdf_random_guesses; 
% % pdfHC = pdfHC./trapz(x, pdfHC);
% 
% % LC pdf
% A = reshape(pLC, size(pLC, 1), 1, size(pLC, 2));
% p_X_LC = p_X.*A;
% pdfLC = sum(p_X_LC, 3);
% pdfLC = pdfLC./trapz(x, pdfLC, 2);
% pdfLC_s = (1 - guessRate)*pdfLC + guessRate*pdf_random_guesses; 
% % pdfLC = pdfLC./trapz(x, pdfLC);
% 
% %% Calculate std dev
% dx = x(2) - x(1); % Assuming uniform
% % warning("Make sure error bins are uniformly spaced")
% mu = 0; %sum( x.*pdf*dx );
% sigma_s = sqrt( sum( ((x - mu).^2).*pdf*dx, 2 ) );
% 
% % For HC
% mu = 0; %sum( x.*pdfHC*dx );
% sigmaHC_s = sqrt( sum( ((x - mu).^2).*pdfHC*dx, 2 ) );
% 
% % For LC
% mu = 0; %sum( x.*pdfLC*dx );
% sigmaLC_s = sqrt( sum( ((x - mu).^2).*pdfLC*dx, 2 ) );
% 
% %% MAD from PDF
% % Calculate metric only after optimization is complete
% if ~optimizationFlag
%     
%     dx = x(2)-x(1);
%     F = cumsum(pdf) * dx;
%     if ~isnan(F(end))
%         F = F / F(end);   % normalize
% %         [Funiq, idx] = unique(F);
% %         xuniq = x(idx);
%         median_val = 0; %interp1(Funiq, xuniq, 0.5);
%         mad_fun = @(d) (interp1(x, F, median_val + d) - interp1(x, F, median_val - d)) - 0.5;
%         
%         d0 = (x(end) - x(1)) / 4; % initial guess for d
%         MAD = fzero(mad_fun, d0);
%         mad_m_s = MAD;
%     else
%         mad_m_s = nan;
%     end
%     
%     % HC
%     dx = x(2)-x(1);
%     F = cumsum(pdfHC) * dx;
%     if ~isnan(F(end))
%         F = F / F(end);   % normalize
% %         [Funiq, idx] = unique(F);
% %         xuniq = x(idx);
%         median_val = 0; %interp1(Funiq, xuniq, 0.5);
%         mad_fun = @(d) (interp1(x, F, median_val + d) - interp1(x, F, median_val - d)) - 0.5;
%         
%         d0 = (x(end) - x(1)) / 4; % initial guess for d
%         MAD = fzero(mad_fun, d0);
%         mad_m_HC_s = MAD;
%     
%     else
%         mad_m_HC_s = nan; % dummy value to avoid nan error
%     end
%     
%     % LC
%     dx = x(2)-x(1);
%     F = cumsum(pdfLC) * dx;
%     if ~isnan(F(end))
%         F = F / F(end);   % normalize
% %         [Funiq, idx] = unique(F);
% %         xuniq = x(idx);
%         median_val = 0; %interp1(Funiq, xuniq, 0.5);
%         mad_fun = @(d) (interp1(x, F, median_val + d) - interp1(x, F, median_val - d)) - 0.5;
%         
%         d0 = (x(end) - x(1)) / 4; % initial guess for d
%         MAD = fzero(mad_fun, d0);
%         mad_m_LC_s = MAD;
%     else
%         mad_m_LC_s = nan;
%     end
% 
% else
%     mad_m_s = 0;
%     mad_m_HC_s = 0;
%     mad_m_LC_s = 0;
% end
% 
% end



% PDF for individual orientations
function [pdf, pdfHC, pdfLC, sigma, sigmaHC, sigmaLC, mad_m, mad_m_HC, mad_m_LC] = getGaussianMixturePDF( ...
    x, sigma_m_stim, bias, pHC, pLC, guessRate, optimizationFlag)
% x = array of perceptual errors

pdf_random_guesses = ones(size(x)) / (numel(x) + 1);
pdf_random_guesses = pdf_random_guesses./trapz(x, pdf_random_guesses); 
pdf_random_guesses = pdf_random_guesses';

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
pdfHC = (1 - guessRate)*pdfHC + guessRate*pdf_random_guesses; 
pdfHC = pdfHC./trapz(x, pdfHC);

% LC pdf
p_X_LC = p_X.*pLC;
pdfLC = sum(p_X_LC, 2);
pdfLC = pdfLC./trapz(x, pdfLC);
pdfLC = (1 - guessRate)*pdfLC + guessRate*pdf_random_guesses; 
pdfLC = pdfLC./trapz(x, pdfLC);

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
% Calculate metric only after optimization is complete
% Compute wrt to zero (assume zero to be mean and the median)
if ~optimizationFlag
    
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
end

end