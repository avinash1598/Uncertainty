% PDF for individual orientations
function [nlls] = getGaussianMixtureNLLs(x, sigma_m_stim, bias, pConf, guessRate, K, confReport)
% x = array of perceptual errors
% P(err|trial) = P(err|trial, confreport)*P(confreport|trial)

% If probability lands outside [-90,90] i.e. period 180, move it back inside by wrapping.
period  = 180;
lls     = zeros(length(x),1);

pConf_norm = pConf./sum(pConf); 

inv2sig2  = 1 ./ (2*sigma_m_stim.^2);
normconst = 1 ./ (sqrt(2*pi) .* sigma_m_stim);

for k = -K:K
    Z = x - bias + k*period;
    p_X = exp(-(Z'.^2).*inv2sig2) .* normconst;
    
    p_X_Conf = p_X.*pConf_norm;
    lls  = lls + sum(p_X_Conf, 2);
end

if confReport == "LC"
    % PDF LC
    pdfRG = 1/period;
    lls = (1 - guessRate)*lls + guessRate*pdfRG;
end

meanProbConf = mean(pConf);
tmp = meanProbConf*lls;
tmp(tmp <= 0) = eps;
nlls = -log( tmp ); % P(conf|ori)*P(err|conf,ori)

end
