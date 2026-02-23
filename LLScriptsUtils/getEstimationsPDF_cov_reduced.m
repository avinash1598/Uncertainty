% Covarying noise

function [retData] = getEstimationsPDF_cov_reduced(rvOriErrs, modelParams, optimizationFlag)

if nargin < 3
    optimizationFlag = false;   % default value
end

sigma_s                   = modelParams.sigma_s;
scale                     = modelParams.scale;
Cc                        = modelParams.Cc;
sigma_meta                = modelParams.sigma_meta;
guessRate                 = modelParams.guessRate;
internalNoiseSamplesCnt   = 1000;

tic
% Internal noise covaries with sensory noise
scaleParam = scale;
shapeParams = sigma_s.^2 ./ scaleParam; 
gammaSamples = zeros(numel(shapeParams), internalNoiseSamplesCnt);

for i = 1:numel(shapeParams)
    shapeParam = shapeParams(i);
    gammaSamples(i, :) = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
        shapeParam, scaleParam);
end

sigma_m = sqrt( gammaSamples ); % Measurement noise

elapsed_time = toc;
disp(['Execution time:::::: ', num2str(elapsed_time), ' seconds']);

% For each value of sigma_m_stim, find the probability of high
% confidence and low confidence. These probabilitites are
% calculated by taking the *inverse* of the log normal distribution
% (which basically is a distribution of confidence variable).
muLogN    = - log((sigma_m.^2) ./ sqrt(sigma_meta.^2 + sigma_m.^2));
sigmaLogN = sqrt(log((sigma_meta.^2)./(sigma_m.^2) + 1));

x1 = repmat(Cc, [1 internalNoiseSamplesCnt]);
cdf_vals  = logncdf(x1, muLogN, sigmaLogN);

pHC_all = 1 - cdf_vals;
pLC_all = cdf_vals;

% Compute PDF for each orientation
[pdf, pdfHC, pdfLC, sigma, sigmaHC, sigmaLC, mad_m, mad_m_HC, mad_m_LC] = getGaussianMixturePDF(rvOriErrs, ...
    sigma_m, pHC_all, pLC_all, guessRate, optimizationFlag);

analyticalPDF    = pdf;
analyticalPDF_HC = pdfHC;
analyticalPDF_LC = pdfLC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Return data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Probability of HC and LC (combined across all orientations)
retData.pHC  = mean(pHC_all);
retData.pLC  = mean(pLC_all);

% Expected sigma - LC and HC - combined across all orientations
% retData.E_sigma_m_HC = sqrt( sum( (sigma_m.^2).*pHC_all, 'all')/sum(pHC_all, 'all') ); % this is right
% retData.E_sigma_m_LC = sqrt( sum( (sigma_m.^2).*pLC_all, 'all')/sum(pLC_all, 'all') );
retData.E_sigma_m_HC = sigmaHC;
retData.E_sigma_m_LC = sigmaLC;

% % Expected sigma for each stimulus orientation - aggreagted by HC and LC:
% Analytical solution

% Expected sigma for all perceptual reports
% retData.E_sigma_m = sqrt( sigma_s.^2 );
retData.E_sigma_m = sigma;

% MAD
retData.mad_m = mad_m;
retData.mad_m_HC = mad_m_HC;
retData.mad_m_LC = mad_m_LC;

% PDF for each orientation
retData.rvOriErrs = rvOriErrs;

% PDF combined across all orientations
retData.analyticalPDF = analyticalPDF;
retData.analyticalPDF_LC = analyticalPDF_LC; % LC
retData.analyticalPDF_HC = analyticalPDF_HC; % HC

end


% PDF for individual orientations
function [pdf, pdfHC, pdfLC, sigma, sigmaHC, sigmaLC, mad_m, mad_m_HC, mad_m_LC] = getGaussianMixturePDF( ...
    x, sigma_m, pHC, pLC, guessRate, optimizationFlag)
% x = array of perceptual errors

pdf_random_guesses = ones(size(x)) / (numel(x) + 1);
pdf_random_guesses = pdf_random_guesses./trapz(x, pdf_random_guesses); 
pdf_random_guesses = pdf_random_guesses';

p_X = exp(- (x').^2 ./ (2*sigma_m.^2 ) ) ./ sqrt(2*pi*sigma_m.^2);
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

mu = 0; %sum( x.*pdf'*dx );
sigma = sqrt( sum( ((x - mu).^2).*pdf'*dx ) );

% For HC
mu = 0; %sum( x.*pdfHC'*dx );
sigmaHC = sqrt( sum( ((x - mu).^2).*pdfHC'*dx ) );

% For LC
mu = 0; %sum( x.*pdfLC'*dx );
sigmaLC = sqrt( sum( ((x - mu).^2).*pdfLC'*dx ) );

%% MAD from PDF
% Calculate metric only after optimization is complete
if ~optimizationFlag
    
    dx = x(2)-x(1);
    F = cumsum(pdf) * dx;
    if ~isnan(F(end))
        F = F / F(end);   % normalize
        %[Funiq, idx] = unique(F);
        %xuniq = x(idx);
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
        %[Funiq, idx] = unique(F);
        %xuniq = x(idx);
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
        %[Funiq, idx] = unique(F);
        %xuniq = x(idx);
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


% keyboard

end