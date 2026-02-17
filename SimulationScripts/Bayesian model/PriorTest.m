
close all
clear all

x = 0:1:180;
mu = 175;
sigma = 10;
likelihood = normpdf(x, mu, sigma);

amp = 0.02;
prior = normpdf(x, 90, 20) + normpdf(x, 0, 20) + normpdf(x, 180, 20);

posterior = likelihood.*prior;
posterior = posterior / trapz(x, posterior);

% for each orientation calculate bias
biases = zeros(1, numel(x));

for i=1:numel(x)
    mu = x(i);
    likelihood = normpdf(x, mu, sigma);
    posterior = likelihood.*prior;
    [val, idx] = max(posterior);
    map_estimate = x(idx);

    bias_ = map_estimate - mu;
    biases(i) = bias_;
end

figure
plot(x, likelihood, DisplayName='likelihood', LineWidth=2)
hold on
plot(x, prior, DisplayName='prior', LineWidth=2)
plot(x, posterior, DisplayName='posterior', LineWidth=2)
plot(x, 0.005*biases, DisplayName="bias scaled")
hold off
legend