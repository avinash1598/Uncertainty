function [retData] = getLLs_SinglyStochastic(trlErrors, trlConfReports, trlStimOris, uniqStimVals, modelParams)

% warning("Deprecated. Use LLScriptsNoBin instead.")

% Precision of orientation error bins
if ~isfield(modelParams, 'oriErrBinWidth')
    modelParams.oriErrBinWidth = 1; %2; % Default value obtained from params check analysis
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
Cc                        = modelParams.Cc;
sigma_meta                = modelParams.sigma_meta;
guessRate                 = modelParams.guessRate;
numOris                   = numel(stimOris);
sigma_m_shape1            = modelParams.sigma_m_shape1;
sigma_m_shape2            = modelParams.sigma_m_shape2;
biasShape                 = modelParams.biasShape;

bias         = biasAmp*sind(biasShape*stimOris);
sigma_s_stim = b + a.*(abs(sind(sigma_m_shape1*stimOris - sigma_m_shape2)));
sigma_m_stim = sigma_s_stim'; % Measurement noise

% For each value of sigma_m_stim, find the probability of high
% confidence and low confidence. These probabilitites are
% calculated by taking the inverse of the log normal distribution
% (which basically is a distribution of confidence variable).
muLogN    = - log((sigma_m_stim.^2) ./ sqrt(sigma_meta.^2 + sigma_m_stim.^2));
sigmaLogN = sqrt(log((sigma_meta.^2)./(sigma_m_stim.^2) + 1));

x1 = repmat(Cc, [numOris 1]);
cdf_vals  = logncdf(x1, muLogN, sigmaLogN);

pHC_all = 1 - cdf_vals;
pLC_all = cdf_vals;

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
        sigma_m_stim(i), bias(i), ...
        pHC_all(i), ...
        guessRate, ...
        K, "HC");

    % LC
    [nllsLC] = getGaussianMixtureNLLs( ...
        fltTrlErrsLC, ...
        sigma_m_stim(i), bias(i), ...
        pLC_all(i), ...
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