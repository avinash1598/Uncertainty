function [retData] = getLLs_ind(trlErrors, trlConfReports, trlStimOris, uniqStimVals, modelParams)

% No of samples to draw for numerical approximation
if ~isfield(modelParams, 'internalNoiseSamplesCnt')
    modelParams.internalNoiseSamplesCnt = 100; % Default value obtained from params check analysis
end

% Precision of orientation error bins
if ~isfield(modelParams, 'oriErrBinWidth')
    modelParams.oriErrBinWidth = 2; % Default value obtained from params check analysis
end

% Controls wrapped gaussian distribution
if ~isfield(modelParams, 'K')
    modelParams.K = 1; % Default value obtained from params check analysis
end

K                         = modelParams.K;
stimOris                  = uniqStimVals;
b                         = modelParams.b;       % Baseline stimulus dependent sensory noise level
a                         = modelParams.a;       % Amplitude for stimulus dependent sensory noise modulation
biasAmp                   = modelParams.biasAmp;
shape                     = modelParams.shape;
scale                     = modelParams.scale;
Cc                        = modelParams.Cc;
sigma_meta                = modelParams.sigma_meta;
guessRate                 = modelParams.guessRate;
internalNoiseSamplesCnt   = modelParams.internalNoiseSamplesCnt; % 100; %1000; %400, 1000;
numOris                   = numel(stimOris);  
sigma_m_shape1            = modelParams.sigma_m_shape1;
sigma_m_shape2            = modelParams.sigma_m_shape2;
biasShape                 = modelParams.biasShape;

bias         = biasAmp*sind(biasShape*stimOris);
sigma_s_stim = b + a.*(abs(sind(sigma_m_shape1*stimOris - sigma_m_shape2)));
% bias = biasAmp*sind(2*stimOris); % As per human subjects
% sigma_s_stim = b + a.*(abs(sind(2*stimOris))); %sigma_s_stim = sigma_s_stim';
% sigma_s_stim = b + a.*(abs(sind(stimOris - 90))); % as per human subjects

% % Internal noise fluctuation (completely independent)
% shapeParam = shape;
% scaleParam = scale;
% gammaSamples = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
%     shapeParam, scaleParam);
% 
% sigma_m_stim = sqrt( sigma_s_stim'.^2 + gammaSamples ); % Measurement noise

% Internal noise fluctuation (completely independent)
% Sample from log normal distribution instead
u = linspace(1/internalNoiseSamplesCnt, 1-1/internalNoiseSamplesCnt, internalNoiseSamplesCnt);

% Find parameters of log normal distribution (i.e. mean and variance of normal distribution I guess)
% mu and sigma here are the parameters (of normal dist) and not necesarily the actual mean
% and variance of the log normal dist
m = shape; v = scale; % desired mean & variance of lognormal distribution
sigma2 = log(1 + v/m^2);  % variance
mu = log(m) - 0.5*sigma2; %mean
lognSamples = logninv(u, mu, sqrt(sigma2));

sigma_m_stim = sqrt( sigma_s_stim'.^2 + lognSamples ); % Measurement noise

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
% mean_cdf_val = mean(cdf_vals(:));
% mean_cdf_val_stim = mean(cdf_vals, 2);

% Compute NLL for each orientation
totalNLL = 0;
trlNLLs  = zeros(1, numel(trlErrors));

for i = 1:numOris
    trlIdx = trlStimOris == stimOris(i);
    trlIdxHC = trlIdx & (trlConfReports == 1);
    trlIdxLC = trlIdx & (trlConfReports == 0);

    fltTrlErrsHC = trlErrors(trlIdxHC); % HC errors 
    fltTrlErrsLC = trlErrors(trlIdxLC); % LC errors

    fltTrlErrsHC = fltTrlErrsHC';
    fltTrlErrsLC = fltTrlErrsLC';

    % HC
    [nllsHC] = getGaussianMixtureNLLs( ...
        fltTrlErrsHC, ...
        sigma_m_stim(i, :), bias(i), ...
        pHC_all(i, :), ...
        guessRate, ...
        K, "HC");

    % LC
    [nllsLC] = getGaussianMixtureNLLs( ...
        fltTrlErrsLC, ...
        sigma_m_stim(i, :), bias(i), ...
        pLC_all(i, :), ...
        guessRate, ...
        K, "LC");

    totalNLL = totalNLL + sum(nllsLC) + sum(nllsHC);
    
    trlNLLs(trlIdxHC) = nllsHC;
    trlNLLs(trlIdxLC) = nllsLC;
end

% return NLL for individual trials as well
retData.totalNLL = totalNLL;
retData.trlNLLs  = trlNLLs;

end
