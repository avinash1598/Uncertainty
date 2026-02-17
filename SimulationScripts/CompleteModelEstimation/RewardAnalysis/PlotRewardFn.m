close all
clear all


%% Analytical reward function
maxTolerableError = 30; %36; % In degrees
x = linspace(-maxTolerableError - 10, maxTolerableError + 10, 2001);

y1 = 18; %18; 15;
x1 = 16; %16.1;
c1 = 0.15; %0.5;
m1 = - c1/y1;
m2 =  ( m1*x1 + c1 ) / (x1 - maxTolerableError); % c2 / maxTolerableError; %m1 - (c1 - c2)/x1;
c2 = - m2*maxTolerableError;

yHC = m1 * sign(x) .* x + c1; yHC(abs(x) > y1) = 1.2*yHC(abs(x) > y1);  yHC(abs(x) > maxTolerableError) = -1.5*c1;
yLC = m2 * sign(x) .* x + c2; yLC(abs(x) > maxTolerableError) = 0; 

x2 = -x1;


scaled_yHC = 100*10*yHC;
scaled_yLC = 100*10*yLC;

yHC_discrete = ceil(scaled_yHC / 10) * 10;
yLC_discrete = ceil(scaled_yLC / 10) * 10;
% yHC_discrete = ceil( (scaled_yHC / 10) * 10);
% yLC_discrete = ceil( (scaled_yLC / 10) * 10);

figure
hold on
plot(x, scaled_yHC, DisplayName="HC", LineWidth=2, Color='green')
plot(x, scaled_yLC, DisplayName="LC", LineWidth=2, Color='red')
plot(x, yHC_discrete, LineWidth=2, Color='green')
plot(x, yLC_discrete, LineWidth=2, Color='red')
xlabel("Perceptual error")
ylabel("Reward")
legend

yline(0, LineWidth=1, LineStyle="--", HandleVisibility="off")
xline(x1, LineWidth=1, LineStyle="--", HandleVisibility="off")
xline(x2, LineWidth=1, LineStyle="--", HandleVisibility="off")

hold off


