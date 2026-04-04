close all
clear; clc;

% Grid
x = linspace(-90,90,2000); % energy should be normalized between -90 and 90

% Gaussian parameters
mu = 0;
sigma = 88;   % make this large to see the difference clearly
% starts to show some affect when it reacher close to 85. Nothing before
% that. K=1 is safe choice option

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


mu = 0;
sigma = 2;

x = -90:0.01:90;
pdf_true = normpdf(x, mu, sigma);

stepSize = [0.001 0.1 0.2];
sampleSizes = [100 1000 10000 100000];

% Decreasing step size too much can overestimte probability in each bin
% This maybe the cause of analytical solution deviating from actual
% solution

figure

for i = 1:length(stepSize)

    N = 100000;

    x = -90:stepSize(i):90;
    pdf_true = normpdf(x, mu, sigma);
    
    % generate samples
    data = mu + sigma*randn(N,1);

    subplot(2,2,i)
    
    % histogram
    histogram(data, -90:stepSize(i):90, "Normalization","pdf")
    hold on
    
    % analytical pdf
    % Trust the analytical solution. PDF estimated from data is impacted as
    % binwidth decreases
    plot(x, pdf_true, "r", "LineWidth",2)

    title(["Step size = " num2str(stepSize(i))])
    xlim([-90 90])

end