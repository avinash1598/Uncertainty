scale = 0.5;

figure
hold on

sigma_e_level = 2;
vals = gamrnd(sigma_e_level/scale, scale, [1 100]);

scatter(1 + 0.1*(rand([1 numel(vals)]) - 0.5), vals)

sigma_e_level = 5;
vals = gamrnd(sigma_e_level/scale, scale, [1 100]);

scatter(2 + 0.1*(rand([1 numel(vals)]) - 0.5), vals)

sigma_e_level = 20;
vals = gamrnd(sigma_e_level/scale, scale, [1 100]);

scatter(3 + 0.1*(rand([1 numel(vals)]) - 0.5), vals)


hold off

% 
% x = [-2, -2, -2, -2, -2];
% % mu = 0;
% % y = (x - mu).^2;
% % v = mean(y);
% 
% % Parameters
% nTrials = 1000;
% sigma_ext = 5;
% sigma_int = 3;
% rho = 0.6; % Desired correlation
% 
% % Covariance matrix
% Sigma = [sigma_ext^2, rho * sigma_ext * sigma_int;
%          rho * sigma_ext * sigma_int, sigma_int^2];
% 
% % Generate bivariate normal noise samples
% noise = mvnrnd([0 0], Sigma, nTrials);  % [eta_ext, eta_int]
% eta_ext = noise(:,1);
% eta_int = noise(:,2);
% 
% % Simulate true orientations and perceived estimates
% theta = rand(nTrials,1) * 180;  % true orientations
% theta_hat = theta + eta_ext + eta_int;  % perceived orientation
% 
% % Check empirical correlation
% empirical_rho = corr(eta_ext, eta_int);
% fprintf('Empirical correlation: %.3f\n', empirical_rho);



% x = zeros(1, 100) - 3 + randn(1, 100);
% mu = 0;
% sigma = sqrt(mean((x - mu).^2)); 
% 
% pdf = normpdf(linspace(-10, 10, 100), mu, sigma);
% 
% figure
% plot(linspace(-10, 10, 100), pdf);
% hold on
% histogram(x, Normalization="probability");
% hold off


% % Parameters for two log-normal distributions
% clear all
% close all
% 
% mu1 = 0.1; sigma1 = 0.5;
% mu2 = 2; sigma2 = 0.5;
% 
% % x-axis values (must be > 0 for log-normal)
% x = linspace(0.01, 10, 500);
% 
% % PDF values
% pdf1 = lognpdf(x, mu1, sigma1);
% pdf2 = lognpdf(x, mu2, sigma2);
% 
% % Plot
% figure; hold on;
% plot(x, pdf1, 'LineWidth', 2, 'DisplayName', '\sigma = 0.5');
% plot(x, pdf2, 'LineWidth', 2, 'DisplayName', '\sigma = 1.0');
% xlabel('x'); ylabel('PDF');
% title('Log-normal Distributions');
% legend; grid on;
