close all
clear all

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Utils')
addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Scripts/')

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Data/COR31.mat');  % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Data/COR32.mat');  % Jiaming
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Data/COR33.mat');  % Akash
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/HumanExpDataAnalysis/Data/CORNFB01.mat'); % Yichao

fltData = data.dat( data.dat.session > 0 , :);
% formattedData = formatExpData(data);
f.dat = fltData;
formattedData = formatExpData(f, false, false);
rvOriErr = -90:3:90;

%% Data structures

orientations               = unique(formattedData.theta_true_all);
n_uncertainty_levels       = formattedData.n_uncertainty_levels;
n_orientations             = numel(orientations);
resp_err_all               = formattedData.resp_err_all;
confidence_report_all      = formattedData.confidence_report_all;
theta_true_all             = formattedData.theta_true_all;

%% Question 1: do variances differ across different conditions
A  = resp_err_all; resp_err_by_level = reshape(A, size(A, 1), []); 
A  = confidence_report_all; confidence_report_by_level   = reshape(A, size(A, 1), []);

% Assumptions
% Independent observations - trials - ok
% Approximately symmetric / not extremely skewed - ok
% Reasonable group sizes - ok
[p, ~] = vartestn(resp_err_by_level', 'TestType', 'BrownForsythe', 'Display', 'on');
% [p, stats] = vartestn(resp_err_by_level', 'TestType', 'LeveneQuadratic', 'Display', 'on');
% [p, stats] = vartestn(resp_err_by_level', 'TestType', 'LeveneAbsolute', 'Display', 'on');

% figure(3)
% % Compute median absolute deviations per condition
% med_abs_dev = zeros(6,1);
% for i = 1:6
%     med_abs_dev(i) = median(abs(resp_err_by_level(i, :) - median( resp_err_by_level(i, :) )));
% end
% bar(med_abs_dev)
% xlabel('Condition')
% ylabel('Median Absolute Deviation')
% title('Brown-Forsythe: group variability')

fprintf("\n")
if p < 0.05
    fprintf('1. Variances differ significantly across uncertainty levels (p = %.4f).\n', p);
else
    fprintf('1. No significant variance differences across uncertainty levels (p = %.4f).\n', p);
end
fprintf("\n")

% %% Question2: Which pair of groups differ in SDs
% pairs = nchoosek(1:n_uncertainty_levels,2);
% p_values = zeros(size(pairs,1),1);
% 
% for k = 1:size(pairs,1)
%     grp1 = resp_err_by_level(pairs(k,1), :); 
%     grp2 = resp_err_by_level(pairs(k,2), :); 
%     values = [grp1(:); grp2(:)];
%     group = [ones(numel(grp1),1); 2*ones(numel(grp2),1)];
%     p_values(k) = vartestn(values, group, 'TestType','BrownForsythe','Display','off');
% end
% 
% % Apply Bonferroni correction
% p_corrected = p_values * size(pairs,1); % Dividing p-values from size would increase the chances of false positives
% p_corrected(p_corrected>1) = 1;
% 
% % Display significant pairs
% for k = 1:size(pairs,1)
%     if p_corrected(k) < 0.05
%         fprintf('Variance differs significantly between condition %d and %d\n', ...
%                  pairs(k,1), pairs(k,2));
%     end
% end


%% Question 2: Do SDs differ by certainty report

results = table('Size',[n_uncertainty_levels 3], ...
                'VariableTypes', {'double','double','string'}, ...
                'VariableNames', {'uncertaintyLevel','pValue','significant'});

for i = 1:n_uncertainty_levels

    e_ = resp_err_by_level(i, :);
    c_ = confidence_report_by_level(i, :);

    values = e_';
    group  = c_';

    p_val = vartestn(values, group, ...
                     'TestType','BrownForsythe', 'Display','off');

    corrected_p = p_val * n_uncertainty_levels;  % Bonferroni

    % Store results
    results.uncertaintyLevel(i) = i;
    results.pValue(i) = corrected_p;

    if corrected_p < 0.05
        results.significant(i) = "YES";
    else
        results.significant(i) = "NO";
    end
end

fprintf("\n")
fprintf('2. Do SDs differ by certainty report?\n');

fprintf("\n")
disp(results)
fprintf('\n')

%% Question 3: Is difference in SDs between low certainty and high certainty because of internal noise?
% Idea of the test: Same stimulus family - does error differes between high
% and low certainty?

errDiff_LCHC = nan + zeros(1, n_uncertainty_levels*n_orientations);
errLCs = nan + zeros(1, n_uncertainty_levels*n_orientations);
errHCs = nan + zeros(1, n_uncertainty_levels*n_orientations);
fltOri = nan + zeros(1, n_uncertainty_levels*n_orientations);

for i=1:n_uncertainty_levels
    for j=1:n_orientations

        errs_         = resp_err_all(i, j, :);
        conf_reports_ = confidence_report_all(i, j, :);
        
        errHC_ = errs_(conf_reports_ == 1);
        errLC_ = errs_(conf_reports_ == 0);
        
        % fprintf('%d, %d \n', numel(errHC_), numel(errLC_))
        
        if numel(errHC_) > 4 && numel(errLC_) > 4 % even with 5 the test is significant

            % For each subsample equal number of values and compute the
            % average error over many permutation
            
            if numel(errHC_) < numel(errLC_)
                y = mean( abs(errHC_) );

                x = 0;
                nitr = 100;
                for itr=1:nitr
                    x_samples = randsample(errLC_(:), numel(errHC_));
                    x = x + mean( abs(x_samples) ) ;
                end
                x = x / nitr;

            elseif numel(errHC_) > numel(errLC_)
                 x = mean( abs(errLC_) );

                 y = 0;
                 nitr = 100;
                 for itr=1:nitr
                    y_samples = randsample(errHC_(:), numel(errLC_));
                    y = y + mean( abs(y_samples) ) ;
                 end
                 y = y / nitr;
            else

                x = mean( abs(errLC_) );
                y = mean( abs(errHC_) );
            end
            
            % x = mean( abs(errLC_) );
            % y = mean( abs(errHC_) );

            errDiff = x - y;
            
            errDiff_LCHC( (i - 1)*n_uncertainty_levels + j ) = errDiff;
            errLCs( (i - 1)*n_uncertainty_levels + j )       = x;
            errHCs( (i - 1)*n_uncertainty_levels + j )       = y;

            fltOri( (i - 1)*n_uncertainty_levels + j )       = orientations(j);

        else
            errDiff = nan;
        end
        
    end
end

[h,p,ci,stats] = ttest(errDiff_LCHC, 0);

fprintf("\n")
if p < 0.05
    fprintf('3. Hint of significant contribution from internal noise (p = %.4f).\n', p);
else
    fprintf('3. No hint of contribution from internal noise  (p = %.4f).\n', p);
end
fprintf("\n")

figure
subplot(2, 2, 1)
hold on
histogram(errDiff_LCHC, 'BinEdges',-40:5:30)
xline(0, 'LineStyle',"--")
xlabel("Abs error diff bw LC and HC")
ylabel("Count")
ylim([0 15])

subplot(2, 2, 2)
hold on
histogram(errLCs, DisplayName="LC") %'BinEdges',0:5:40, 
histogram(errHCs, DisplayName="HC") % 'BinEdges',0:5:40,
xline(0, 'LineStyle',"--")
xlabel("Abs error")
ylabel("Count")
legend


%% Marginal distribution: Question 4: do variances differ across different contrast level?
resp_err_by_c_level          = formattedData.marginal_contrast_resp_err_all;

[p, ~] = vartestn(resp_err_by_c_level', 'TestType', 'BrownForsythe', 'Display', 'on');

fprintf("\n")
if p < 0.05
    fprintf('4. Variances differ significantly across contrast levels (p = %.4f).\n', p);
else
    fprintf('4. No significant variance differences across contrast levels (p = %.4f).\n', p);
end
fprintf("\n")

%% Marginal distribution: Question 5: do variances differ across different spread amount?
resp_err_by_s_level          = formattedData.marginal_spread_resp_err_all;

[p, stats] = vartestn(resp_err_by_s_level', 'TestType', 'BrownForsythe', 'Display', 'on');

fprintf("\n")
if p < 0.05
    fprintf('5. Variances differ significantly across spread amount (p = %.4f).\n', p);
else
    fprintf('5. No significant variance differences across spread amount (p = %.4f).\n', p);
end
fprintf("\n")


%% Marginal distribution: Question 6: do variances differ across different duration?
resp_err_by_d_level          = formattedData.marginal_dur_resp_err_all;

[p, stats] = vartestn(resp_err_by_d_level', 'TestType', 'BrownForsythe', 'Display', 'on');

fprintf("\n")
if p < 0.05
    fprintf('6. Variances differ significantly across stim duration (p = %.4f).\n', p);
else
    fprintf('6. No significant variance differences across stim duration (p = %.4f).\n', p);
end
fprintf("\n")

%% Marginal distribution: Question 7: Do SDs differ by certainty report for contrast manipulation? 

c_levels                     = formattedData.marginal_contrast_vals;
resp_err_by_c_level          = formattedData.marginal_contrast_resp_err_all;
confidence_report_by_c_level = formattedData.marginal_contrast_confidence_report_all;

results = table('Size',[numel(c_levels) 3], ...
                'VariableTypes', {'double','double','string'}, ...
                'VariableNames', {'contrast','pValue','significant'});

for i = 1:numel(c_levels)

    e_ = resp_err_by_c_level(i, :);
    c_ = confidence_report_by_c_level(i, :);

    values = e_';
    group  = c_';

    p_val = vartestn(values, group, ...
                     'TestType','BrownForsythe', 'Display','off');

    corrected_p = p_val * numel(c_levels);  % Bonferroni

    % Store results
    results.contrast(i) = c_levels(i);
    results.pValue(i) = corrected_p;

    if corrected_p < 0.05
        results.significant(i) = "YES";
    else
        results.significant(i) = "NO";
    end
end

fprintf("\n")
fprintf('7. Do SDs differ by certainty report for contrast manipulation? \n');

fprintf("\n")
disp(results)
fprintf('\n')

%% Marginal distribution: Question 8: Do SDs differ by certainty report for spread manipulation? 

s_levels                     = formattedData.marginal_spread_vals2;
resp_err_by_s_level          = formattedData.marginal_spread_resp_err_all2;
confidence_report_by_s_level = formattedData.marginal_spread_confidence_report_all2;

results = table('Size',[numel(s_levels) 3], ...
                'VariableTypes', {'double','double','string'}, ...
                'VariableNames', {'spread','pValue','significant'});

for i = 1:numel(s_levels)

    e_ = resp_err_by_s_level(i, :);
    c_ = confidence_report_by_s_level(i, :);

    values = e_';
    group  = c_';

    p_val = vartestn(values, group, ...
                     'TestType','BrownForsythe', 'Display','off');

    corrected_p = p_val * numel(s_levels);  % Bonferroni

    % Store results
    results.spread(i) = s_levels(i);
    results.pValue(i) = corrected_p;

    if corrected_p < 0.05
        results.significant(i) = "YES";
    else
        results.significant(i) = "NO";
    end
end

fprintf("\n")
fprintf('8. Do SDs differ by certainty report for spread manipulation? \n');

fprintf("\n")
disp(results)
fprintf('\n')

%% Marginal distribution: Question 9: Do SDs differ by certainty report for duration manipulation? 

d_levels                     = formattedData.marginal_dur_vals;
resp_err_by_d_level          = formattedData.marginal_dur_resp_err_all;
confidence_report_by_d_level = formattedData.marginal_dur_confidence_report_all;

results = table('Size',[numel(d_levels) 3], ...
                'VariableTypes', {'double','double','string'}, ...
                'VariableNames', {'duration','pValue','significant'});

for i = 1:numel(d_levels)

    e_ = resp_err_by_d_level(i, :);
    c_ = confidence_report_by_d_level(i, :);

    values = e_';
    group  = c_';

    p_val = vartestn(values, group, ...
                     'TestType','BrownForsythe', 'Display','off');

    corrected_p = p_val * numel(d_levels);  % Bonferroni

    % Store results
    results.duration(i) = d_levels(i);
    results.pValue(i) = corrected_p;

    if corrected_p < 0.05
        results.significant(i) = "YES";
    else
        results.significant(i) = "NO";
    end
end

fprintf("\n")
fprintf('9. Do SDs differ by certainty report for duration manipulation? \n');

fprintf("\n")
disp(results)
fprintf('\n')


%% Question 10: do perceptual error show signature of multiplicative effect?

% Bootstrap test done here might not be right way to do this. x-axis should
% probablt be the sensory noise. But sensory noise can only be obtained
% after model fitting! So maybe model fitting is probably the right way to
% differentiate between alternative hypothesis. 

figure(16)

% x = mean(resp_err_by_level, 2);
x = median(resp_err_by_level, 2);
% y = std(resp_err_by_level, 0, 2);
y = mad(resp_err_by_level, 1, 2);

HC_idx = confidence_report_by_level == 1;
LC_idx = confidence_report_by_level == 0;

resp_HC = resp_err_by_level;
resp_HC(~HC_idx) = NaN;

resp_LC = resp_err_by_level;
resp_LC(~LC_idx) = NaN;

x_HC = mean(resp_HC, 2, 'omitnan');
% y_HC = std(resp_HC, 0, 2, 'omitnan');
y_HC = mad(resp_HC, 1, 2);
% y_HC = y_HC.^2;

x_LC = mean(resp_LC, 2, 'omitnan');
% y_LC = std(resp_LC, 0, 2, 'omitnan');
y_LC = mad(resp_LC, 1, 2);
% y_LC = y_LC.^2;

subplot(3, 3, 1)
errorbar(1:n_uncertainty_levels, ...
    x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6);

xlabel("Uncertainty level")
ylabel("Std(data)")

subplot(3, 3, 2)

% Behavioral variability
plot(1:n_uncertainty_levels, y, LineWidth=2);
xlabel("Uncertainty level")
ylabel("Std(data)")
xticks(1:n_uncertainty_levels);
xticklabels(1:n_uncertainty_levels);

subplot(3, 3, 3)

x_ = y; %1:n_uncertainty_levels;
y_ = y_HC;
slopeHC = polyfit(x_, y_, 1);
xfit_HC = linspace(min(x_), max(x_), 100);
yfit_HC = polyval(slopeHC, xfit_HC);

x_ = y; %1:n_uncertainty_levels;
y_ = y_LC;
slopeLC = polyfit(x_, y_, 1);
xfit_LC = linspace(min(x_), max(x_), 100);
yfit_LC = polyval(slopeLC, xfit_LC);


% Behavioral variability
hold on
scatter(y, y_HC,  LineWidth=2, DisplayName="HC");
plot(xfit_HC, yfit_HC, 'LineWidth', 2);
scatter(y, y_LC,  LineWidth=2, DisplayName="LC");
plot(xfit_LC, yfit_LC, 'LineWidth', 2);
legend
% xticks(1:n_uncertainty_levels);
% xticklabels(1:n_uncertainty_levels);
xlabel("Uncertainty level")
ylabel("Std(data)")
hold off


% should this be varaince or std
uncertainty_level_all = repmat((1:n_uncertainty_levels)', 1, size(resp_err_by_level, 2));
error                 = abs(resp_err_by_level(:)) + 1e-12;
uncertaintyLevel      = uncertainty_level_all(:);
confidence            = confidence_report_by_level(:);
T                     = table(error, uncertaintyLevel, confidence);
T.confidence          = categorical(T.confidence);  % HC vs LC

% % lm = fitlm(T, 'error ~ uncertaintyLevel * confidence');
% 
% % This analysis might be okay if i correct for bias
% glme = fitglme(T, ...
%     'error ~ uncertaintyLevel * confidence + (1|uncertaintyLevel)', ...
%     'Distribution','Gamma', ...
%     'Link','log');
% 
% disp(glme)


err_HC_1 = resp_err_by_level(1, confidence_report_by_level(1, :) == 1);
err_HC_2 = resp_err_by_level(6, confidence_report_by_level(6, :) == 1);
err_LC_1 = resp_err_by_level(1, confidence_report_by_level(1, :) == 0);
err_LC_2 = resp_err_by_level(6, confidence_report_by_level(6, :) == 0);

% d1 = std(err_LC_1) - std(err_HC_1);
% d2 = std(err_LC_2) - std(err_HC_2);
d1 = mad(err_LC_1, 1) - mad(err_HC_1, 1);
d2 = mad(err_LC_2, 1) - mad(err_HC_2, 1);
D_obs = d2 - d1;

nBoot = 10000;
D_boot = zeros(nBoot,1);

% Do this with simulated data

for b = 1:nBoot
%     sHC1 = std(err_HC_1(randi(numel(err_HC_1), [numel(err_HC_1), 1])));
%     sLC1 = std(err_LC_1(randi(numel(err_LC_1), [numel(err_LC_1), 1])));
%     sHC2 = std(err_HC_2(randi(numel(err_HC_2), [numel(err_HC_2), 1])));
%     sLC2 = std(err_LC_2(randi(numel(err_LC_2), [numel(err_LC_2), 1])));
    
    sHC1 = mad( err_HC_1(randi(numel(err_HC_1), [numel(err_HC_1), 1])), 1);
    sLC1 = mad( err_LC_1(randi(numel(err_LC_1), [numel(err_LC_1), 1])), 1);
    sHC2 = mad( err_HC_2(randi(numel(err_HC_2), [numel(err_HC_2), 1])), 1);
    sLC2 = mad( err_LC_2(randi(numel(err_LC_2), [numel(err_LC_2), 1])), 1);
    
    D_boot(b) = (sLC2 - sHC2) - (sLC1 - sHC1);
end

p = mean(D_boot >= D_obs);

% 95% percentile CI
ci = prctile(D_boot, [2.5 97.5]);

% one-sided p estimate (probability D_boot <= 0)
p_one_sided = mean(D_boot <= 0); % small means evidence D>0

fprintf('D_obs=%.4f, CI=[%.4f, %.4f], p_one_sided (D>0) ~= %.4f\n', D_obs, ci(1), ci(2), p_one_sided);

figure
hold on
histogram(D_boot)
xline(D_obs, LineStyle="--")
xline(ci(1), LineStyle="-")
xline(ci(2), LineStyle="-")
hold off

retData = doBootstrapTest(abs(err_HC_1), abs(err_LC_1), abs(err_HC_2), abs(err_LC_2), true);
retData

%% Marginals

c_levels                     = formattedData.marginal_contrast_vals;
resp_err_by_c_level          = formattedData.marginal_contrast_resp_err_all;
confidence_report_by_c_level = formattedData.marginal_contrast_confidence_report_all;

err_HC_1 = resp_err_by_c_level(1, confidence_report_by_c_level(1, :) == 1);
err_HC_2 = resp_err_by_c_level(2, confidence_report_by_c_level(2, :) == 1);
err_LC_1 = resp_err_by_c_level(1, confidence_report_by_c_level(1, :) == 0);
err_LC_2 = resp_err_by_c_level(2, confidence_report_by_c_level(2, :) == 0);

retData = doBootstrapTest(abs(err_HC_1), abs(err_LC_1), abs(err_HC_2), abs(err_LC_2), true);
retData


s_levels                     = formattedData.marginal_spread_vals2;
resp_err_by_s_level          = formattedData.marginal_spread_resp_err_all2;
confidence_report_by_s_level = formattedData.marginal_spread_confidence_report_all2;

err_HC_1 = resp_err_by_s_level(1, confidence_report_by_s_level(1, :) == 1);
err_HC_2 = resp_err_by_s_level(2, confidence_report_by_s_level(2, :) == 1);
err_LC_1 = resp_err_by_s_level(1, confidence_report_by_s_level(1, :) == 0);
err_LC_2 = resp_err_by_s_level(2, confidence_report_by_s_level(2, :) == 0);

retData = doBootstrapTest(abs(err_HC_1), abs(err_LC_1), abs(err_HC_2), abs(err_LC_2), true);
retData


d_levels                     = formattedData.marginal_dur_vals;
resp_err_by_d_level          = formattedData.marginal_dur_resp_err_all;
confidence_report_by_d_level = formattedData.marginal_dur_confidence_report_all;

err_HC_1 = resp_err_by_d_level(1, confidence_report_by_d_level(1, :) == 1);
err_HC_2 = resp_err_by_d_level(2, confidence_report_by_d_level(2, :) == 1);
err_LC_1 = resp_err_by_d_level(1, confidence_report_by_d_level(1, :) == 0);
err_LC_2 = resp_err_by_d_level(2, confidence_report_by_d_level(2, :) == 0);

retData = doBootstrapTest(abs(err_HC_1), abs(err_LC_1), abs(err_HC_2), abs(err_LC_2), true);
retData
