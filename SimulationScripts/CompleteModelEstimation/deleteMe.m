close all
clear; clc;

% Grid
x = linspace(-90,90,2000); % energy should be normalized between -90 and 90

% Gaussian parameters
mu = 0;
sigma = 40;   % make this large to see the difference clearly

% ---- Standard Gaussian (infinite support) ----
gauss = exp(-(x-mu).^2 ./ (2*sigma^2)) ./ sqrt(2*pi*sigma^2);

% ---- Truncated Gaussian (your method) ----
trunc = gauss ./ trapz(x,gauss);

% ---- Wrapped Gaussian ----
wrap = zeros(size(x));

% number of wraps to approximate infinite sum
K = 2;
temp = zeros(2*K+1, numel(x));

for k = -K:K
    temp(k+K+1, :) = exp(-(x - mu + 180*k).^2 ./ (2*sigma^2)) ./ sqrt(2*pi*sigma^2);
    wrap = wrap + exp(-(x - mu + 180*k).^2 ./ (2*sigma^2)) ./ sqrt(2*pi*sigma^2);
end

% normalize wrapped distribution
wrap = wrap ./ trapz(x,wrap);

% ---- Plot ----
figure
plot(x,gauss,'k--','LineWidth',2); hold on
plot(x,trunc,'r','LineWidth',2)
plot(x,wrap,'b','LineWidth',2)

xlabel('Orientation error (deg)')
ylabel('Probability density')

legend('Original Gaussian','Truncated Gaussian','Wrapped Gaussian')

% xlim([-90 90])
grid on

figure
plot(x, temp, LineWidth=2)


% mu = 0;
% sigma = 10;
% 
% x = -90:0.1:90;
% pdf_true = normpdf(x, mu, sigma);
% 
% sampleSizes = [100 1000 10000 100000];
% 
% figure
% 
% for i = 1:length(sampleSizes)
% 
%     N = sampleSizes(i);
% 
%     % generate samples
%     data = mu + sigma*randn(N,1);
% 
%     subplot(2,2,i)
% 
%     % histogram
%     histogram(data, -90:0.1:90, "Normalization","pdf")
%     hold on
% 
%     % analytical pdf
%     plot(x, pdf_true, "r", "LineWidth",2)
% 
%     title(["N = " num2str(N)])
%     xlim([-90 90])
% 
% end