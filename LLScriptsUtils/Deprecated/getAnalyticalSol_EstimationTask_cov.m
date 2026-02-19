% Covarying noise

function [retData] = getAnalyticalSol_EstimationTask_cov(modelParams)

stimOris                  = modelParams.orientations;
b                         = modelParams.b;       % Baseline stimulus dependent sensory noise level
a                         = modelParams.a;       % Amplitude for stimulus dependent sensory noise modulation
biasAmp                   = modelParams.biasAmp;
scale                     = modelParams.scale;
Cc                        = modelParams.Cc;
sigma_meta                = modelParams.sigma_meta;
internalNoiseSamplesCnt   = 1000;
numOris                   = numel(stimOris);

sigma_s_stim = b + a.*(abs(sind(2*stimOris))); %sigma_s_stim = sigma_s_stim';

% Internal noise covaries with sensory noise
shapeParams = sigma_s_stim;
scaleParam  = scale;

gammaSamples = zeros(numel(shapeParams), internalNoiseSamplesCnt);

for i = 1:numel(shapeParams)
    shapeParam = shapeParams(i);
    scaleParam = scale;
    gammaSamples(i, :) = gaminv(linspace(1/internalNoiseSamplesCnt, 1 - 1/internalNoiseSamplesCnt, internalNoiseSamplesCnt), ...
        shapeParam, scaleParam);
end

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
% I1 = trapz(gammaSamples, y, 2); % Integrate y with non-uniform spacing gammaSamples
% I2 = trapz(stimOris, I1*(1/numOris));
dx = gammaSamples(:,2:end) - gammaSamples(:,1:end-1);                          
avgY = (y(:,1:end-1) + y(:,2:end)) / 2;                  
I1 = sum(dx .* avgY, 2);                                 

% Is I2 correct way to find the expectation of PLC? I2 might be equal to
% mean_cdf_val only when no of thetas is large. Run simulations to see what
% is analytical solution. I would assume mean_cdf_val to be the correct solution.

diff = abs(mean_cdf_val_stim - I1);
assert(all(diff < 0.2), 'Sanity check failed: not all differences < 0.1'); % Make sure the two values are same (just a sanity check for WLLN)

% if sum(diff > 0.2) > 0
%     disp(diff)
%     size(diff)
% end

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
retData.bias = biasAmp*sind(2*stimOris);

end



% Returns
% Expected sigma (aggregate and by orientation)
% PHC and PLC (aggregate and by orientation)
% bias (although this might not be needed in this function)
% TODO: maybe also return distribution for each orientation. bias will be
% used in that case. Also return distribution merged across all the
% orientations


% Expecyted sigma - aggreagted by HC and LC
% y = sigma_m_stim.*pdf;
% retData.E_sigma_m_stim_by_ori = trapz(gammaSamples, y, 2);
% retData.E_sigma_m_stim_by_ori = mean(sigma_m_stim, 2);


% figure
% subplot(2, 3, 1)
% plot(stimOris, sigma_s_stim.^2)
% subplot(2, 3, 2)
% plot(gammaSamples, gamcdf(gammaSamples, shapeParam, scaleParam))
% hold on
% plot(gammaSamples, gampdf(gammaSamples, shapeParam, scaleParam))
% hold off
% subplot(2, 3, 3)
% vals = ( mean(sigma_m_stim.^2, 2) ); % Matches analytical solution (1.99, 2.91)
% plot(stimOris, vals)

% subplot(2, 3, 4)
% plot(stimOris, I1)
% hold on
% plot(stimOris, retData.pLC_byOri) % These two are literally same when no of samples are large
% hold off
% 
% subplot(2, 3, 5)
% plot(stimOris, retData.bias)
