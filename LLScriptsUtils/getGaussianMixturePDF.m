function p=getGaussianMixturePDF(x, sigma_s, sigma_si, varGain, bias)

%%%
% x = array of perceptual errors
%%%

sampleCnt = 500;
shapeParam = 1./varGain;
scaleParam = varGain;
gammaSamples = gaminv(linspace(1/sampleCnt, 1 - 1/sampleCnt, sampleCnt), ...
    shapeParam, scaleParam);
sigma_m = sqrt( sigma_s'.^2 + (sigma_si*gammaSamples).^2 );

x = x - bias;
p_X = exp(- (x').^2 ./ (2*sigma_m.^2 ) ) ./ sqrt(2*pi*sigma_m.^2);
p = sum(p_X, 2);
p = p./trapz(x, p);

end


% % Average PDF over all the orientation to get analytical solution
% x = -90:0.5:90;
% analyticalPDF = zeros(1, numel(x));
% for i=1:n_theta
%     pdf_ = getGaussianMixturePDF(x, sigma_s_stim(i), sigma_si, varGain, bias(i));
%     analyticalPDF = analyticalPDF + pdf_;
% end
% analyticalPDF = analyticalPDF/n_theta;
% analyticalPDF = analyticalPDF / trapz(x, analyticalPDF);
