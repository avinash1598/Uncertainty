close all
clear all


sigma_theta = 10;

sigma_reward = 5; % start with sigma reward function
eps = linspace(-40, 40, 1000); % error values

R_HC = normpdf(eps, 0, sigma_reward);        % Reward function
R_LC = 10*ones(1, length(eps)) / length(eps);
P = normpdf(eps, 0, sigma_theta);        % Error distribution

% Compute convolution
conv_RP_HC = conv(R_HC, P, 'same');
conv_RP_LC = conv(R_LC, P, 'same');


figure

subplot(2, 2, 1)
hold on
plot(eps, R_LC, DisplayName="Reward LC")
plot(eps, R_HC, DisplayName="Reward HC")
plot(eps, P, DisplayName="Error distribution")
xlabel("Error")
ylabel("Reward amount")
legend()
hold off

% Plot
subplot(2, 2, 2)
hold on
plot(eps, conv_RP_HC, 'DisplayName', "HC");
plot(eps, conv_RP_LC, 'DisplayName', "LC");

xlabel('Reported error');
ylabel('Expected reward');
% title('Expected reward with error');
legend()
xlim([0, 40])
hold off


% Expected reward at different orientation error
sigma_theta = [5, 10, 15, 20, 25];

sigma_reward = 5; % start with sigma reward function
eps = linspace(-40, 40, 1000); % error values

R_HC = normpdf(eps, 0, sigma_reward);          % Reward function
R_LC = 10*ones(1, length(eps)) / length(eps);

subplot(2, 2, 3)
hold on

for idx=1:length(sigma_theta)
    % Error distribution
    P = normpdf(eps, 0, sigma_theta(idx)); 

    % Compute convolution
    conv_RP_HC = conv(R_HC, P, 'same');
    conv_RP_LC = conv(R_LC, P, 'same');

    plot(eps, conv_RP_HC, 'DisplayName', sprintf('\\sigma_p = %.1f', sigma_theta(idx)));
end

xlabel('Reported error');
ylabel('Expected reward');
legend()
ylim([0, 1])
xlim([0, 40])
hold off


subplot(2, 2, 4)
hold on

for idx=1:length(sigma_theta)
    % Error distribution
    P = normpdf(eps, 0, sigma_theta(idx)); 
    
    % Compute convolution
    conv_RP_HC = conv(R_HC, P, 'same');
    conv_RP_LC = conv(R_LC, P, 'same');

    plot(eps, conv_RP_LC, 'DisplayName', sprintf('\\sigma_p = %.1f', sigma_theta(idx)));
end

xlabel('Reported error');
ylabel('Expected reward');
legend()
xlim([0, 40])
ylim([0, 1])
hold off