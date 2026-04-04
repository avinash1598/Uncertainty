close all
clear all

% desired mean/variance of lognormal
m = 10; v = 25;

% convert to log-space params
sigma2 = log(1 + v/m^2);
mu = log(m) - 0.5*sigma2;

% sample
N = 1e5;
x = lognrnd(mu, sqrt(sigma2), [N,1]);

% plot + pdf
figure
histogram(x, 100, 'Normalization','pdf'); hold on
xx = linspace(min(x), max(x), 200);
plot(xx, lognpdf(xx, mu, sqrt(sigma2)), 'r','LineWidth',2)