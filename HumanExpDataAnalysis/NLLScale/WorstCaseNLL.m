restoredefaultpath
% close all
clear all

set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',1, ...
    'defaultAxesBox','off', ...
    'defaultAxesTickDir','out', ...
    'defaultLineLineWidth',1.5, ...
    'defaultAxesLabelFontSizeMultiplier',1.1, ...
    'defaultAxesTitleFontWeight','normal', ...
    'defaultLegendBox','off');


addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Utils/')

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR31.mat'); % Tien
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/COR33.mat'); % Akash
data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB01.mat'); % Yichao
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB02.mat');   % Jonathan
% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/ProcessModel/HumanExpDataAnalysis/Data/CORNFB03.mat'); % David

% data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR32.mat'); % Jiaming

fltData = data.dat( data.dat.session > 0, :);
f.dat = fltData;
formattedData = formatExpData(f, false, true);
rvOriErr = -90:2:90; % Why does NLL decreases as bin becomes very fine?

orientations               = unique(formattedData.theta_true_all);
resp_err_all               = formattedData.resp_err_all;
confidence_report_all      = formattedData.confidence_report_all;
resp_err_all_flat          = formattedData.resp_err_all(:);
confidence_report_all_flat = formattedData.confidence_report_all(:);

n_uncertainty_levels = formattedData.n_uncertainty_levels;
A = resp_err_all; resp_err_by_level = reshape(A, size(A, 1), []); % is this correct
A = confidence_report_all; confidence_report_by_level = reshape(A, size(A, 1), []);


%% Both error distribution and confidence report random
N = numel(resp_err_all);
error_range = 90 - (-90);    % 180 degrees
P_e = 1 / error_range;
P_c = 0.5;
P_stim = 1/numel(orientations);
NLL_uniform_by_stim = -N * log(P_e * P_c ); % * P_stim
NLL_uniform = -N * log(P_e * P_c);
% disp(NLL_uniform);

fprintf('NLL: %.4f \n', NLL_uniform);
fprintf('NLL by ori: %.4f \n', NLL_uniform_by_stim);

%% Error distribution not random but confidence report random

% Fit PDFs to shuffled distribution
% Calculate NLL on actual data

% Shuffle confidence label
confidence_report_by_level_shuffled = confidence_report_by_level( ...
    randperm(numel(confidence_report_by_level)));
confidence_report_by_level_shuffled = reshape( ...
    confidence_report_by_level_shuffled, ...
    size(confidence_report_by_level));

% resp_err_by_level_shuffled = resp_err_by_level( ...
%     randperm(numel(resp_err_by_level)));
% resp_err_by_level_shuffled = reshape( ...
%     resp_err_by_level_shuffled, ...
%     size(resp_err_by_level));

totalNLL = 0;

figure

confLabels = ["LC", "HC"];
for confVal = [0, 1]
    for i=1:n_uncertainty_levels

        subplot(3, n_uncertainty_levels, i + n_uncertainty_levels*confVal)
        hold on

        % oriErrByLevel_shuffled    = resp_err_by_level_shuffled(i, :);
        oriErrByLevel    = resp_err_by_level(i, :);
        reportedConf     = confidence_report_by_level(i, :);
        reportedConf_shuffled = confidence_report_by_level_shuffled(i, :);
        grpOriErr          = oriErrByLevel(reportedConf == confVal);
        grpOriErr_shuffled = oriErrByLevel(reportedConf_shuffled == confVal);
        % grpOriErr_shuffled = oriErrByLevel_shuffled(reportedConf_shuffled == confVal);
        
        pConfReport = sum(reportedConf == confVal) / numel(reportedConf);
        pConfReport_shuffled = sum(reportedConf_shuffled == confVal) / numel(reportedConf_shuffled);
        % disp(pConfReport_shuffled)

        % Fit PDF on shuffled data
        [~, xi] = ksdensity(grpOriErr_shuffled, 'Function', 'pdf', 'Support', [-90 90]);
        
        % Hack to avoid NaN issue because of shuffling
        dx = xi(2) - xi(1);
        edges = -90:dx:90;
        if edges(end) < 90
            edges = [edges, edges(end)+dx]; % extend so ≥90
        end
        
        % compute KDE again
        [f, xi] = ksdensity(grpOriErr_shuffled, edges, 'Function', 'pdf', 'Support', [-90 90]); % 'Support', [-90 90]
        
        % disp(trapz(xi, f))
        %f = f ./ trapz(xi, f);
        
        % Calculate likelihood on actual unshuffled data
        likelihood = calculateLikelihood(f, xi, grpOriErr);
        %assert(sum(isnan(likelihood)) == 0)
        % likelihood = likelihood.*pConfReport;        %
        likelihood = likelihood.*pConfReport_shuffled; % Model here is fit on shuffled data so pConf should come from shuffled data
        nll_ = -log(likelihood);
        totalNLL = totalNLL + sum(nll_);
        
        % Plot
        histogram(grpOriErr, rvOriErr, Normalization="pdf", DisplayName=confLabels(confVal + 1));
        hold on
        plot(xi, f, DisplayName="PDF", LineWidth=1.5)
        hold off
        
        xlabel("Orientation (deg)")
        ylabel("count")
        title(sprintf('C: %.2f D: %.2f T: %.2f ', ...
            formattedData.uncertaintyVals(i, 1), ...
            formattedData.uncertaintyVals(i, 2), ...
            formattedData.uncertaintyVals(i, 3)))
        legend

    end
end

fprintf('NLL (shuffled conf label): %.4f \n', totalNLL);

% comaprison fo reduced vs full might not be a good idea
% 70982 - 7541
% Each bin - prob

%% Uncertainty level and confidence report random

% Fit PDFs to shuffled distribution
% Calculate NLL on actual data

% Shuffle confidence label
confidence_report_by_level_shuffled = confidence_report_by_level( ...
    randperm(numel(confidence_report_by_level)));
confidence_report_by_level_shuffled = reshape( ...
    confidence_report_by_level_shuffled, ...
    size(confidence_report_by_level));

resp_err_by_level_shuffled = resp_err_by_level( ...
    randperm(numel(resp_err_by_level)));
resp_err_by_level_shuffled = reshape( ...
    resp_err_by_level_shuffled, ...
    size(resp_err_by_level));

totalNLL = 0;

figure

confLabels = ["LC", "HC"];
for confVal = [0, 1]
    for i=1:n_uncertainty_levels

        subplot(3, n_uncertainty_levels, i + n_uncertainty_levels*confVal)
        hold on

        oriErrByLevel_shuffled    = resp_err_by_level_shuffled(i, :);
        oriErrByLevel    = resp_err_by_level(i, :);
        reportedConf = confidence_report_by_level(i, :);
        reportedConf_shuffled = confidence_report_by_level_shuffled(i, :);
        grpOriErr    = oriErrByLevel(reportedConf == confVal);
        % grpOriErr_shuffled = oriErrByLevel(reportedConf_shuffled == confVal);
        grpOriErr_shuffled = oriErrByLevel_shuffled(reportedConf_shuffled == confVal);
        
        pConfReport = sum(reportedConf == confVal) / numel(reportedConf);
        pConfReport_shuffled = sum(reportedConf_shuffled == confVal) / numel(reportedConf_shuffled);
        % disp(pConfReport_shuffled)

        % Fit PDF on shuffled data
        [~, xi] = ksdensity(grpOriErr_shuffled, 'Function', 'pdf', 'Support', [-90 90]);
        
        % Hack to avoid NaN issue because of shuffling
        dx = xi(2) - xi(1);
        edges = -90:dx:90;
        if edges(end) < 90
            edges = [edges, edges(end)+dx]; % extend so ≥90
        end
        
        % compute KDE again
        [f, xi] = ksdensity(grpOriErr_shuffled, edges, 'Function', 'pdf', 'Support', [-90 90]); % 'Support', [-90 90]
        
        % disp(trapz(xi, f))
        %f = f ./ trapz(xi, f);
        
        % Calculate likelihood on actual unshuffled data
        likelihood = calculateLikelihood(f, xi, grpOriErr);
        %assert(sum(isnan(likelihood)) == 0)
        % likelihood = likelihood.*pConfReport;        %
        likelihood = likelihood.*pConfReport_shuffled; % Model here is fit on shuffled data so pConf should come from shuffled data
        nll_ = -log(likelihood);
        totalNLL = totalNLL + sum(nll_);
        
        % Plot
        histogram(grpOriErr, rvOriErr, Normalization="pdf", DisplayName=confLabels(confVal + 1));
        hold on
        plot(xi, f, DisplayName="PDF", LineWidth=1.5)
        hold off
        
        xlabel("Orientation (deg)")
        ylabel("count")
        title(sprintf('C: %.2f D: %.2f T: %.2f ', ...
            formattedData.uncertaintyVals(i, 1), ...
            formattedData.uncertaintyVals(i, 2), ...
            formattedData.uncertaintyVals(i, 3)))
        legend

    end
end

fprintf('NLL (shuffled conf label and uncertainty level): %.4f \n', totalNLL);

%% By ori
confidence_report_all_shuffled = confidence_report_all( ...
    randperm(numel(confidence_report_all)));
confidence_report_all_shuffled = reshape( ...
    confidence_report_all_shuffled, ...
    size(confidence_report_all));

resp_err_all_shuffled = resp_err_all( ...
    randperm(numel(resp_err_all)));
resp_err_all_shuffled = reshape( ...
    resp_err_all_shuffled, ...
    size(resp_err_all));

% % Sanity test: equivalent to best case
% confidence_report_all_shuffled = confidence_report_all;
% resp_err_all_shuffled = resp_err_all;

totalNLL = 0;

for confVal = [0, 1]
    for i=1:n_uncertainty_levels
        for j=1:numel(orientations)

            oriErr_shuffled  = squeeze( resp_err_all_shuffled(i, j, :) );
            oriErr           = squeeze( resp_err_all(i, j, :) );
            reportedConf     = squeeze( confidence_report_all(i, j, :) );
            reportedConf_shuffled = squeeze( confidence_report_all_shuffled(i, j, :) );
            grpOriErr        = oriErr(reportedConf == confVal);
            %grpOriErr_shuffled = oriErr(reportedConf_shuffled == confVal);
            grpOriErr_shuffled = oriErr_shuffled(reportedConf_shuffled == confVal);
            
            pConfReport = sum(reportedConf == confVal) / numel(reportedConf);
            pConfReport_shuffled = sum(reportedConf_shuffled == confVal) / numel(reportedConf_shuffled);
            %pStim = 1/numel(orientations);
            
            % P(err, conf, stim) = P(err|conf,stim)P(conf|stim)P(stim)
            
            if numel(grpOriErr) == 0
                continue
            end
            
            %assert(numel(grpOriErr) > 0);
            
            % Fit PDF on shuffled data
            [~, xi] = ksdensity(grpOriErr_shuffled, 'Function', 'pdf', 'Support', [-90 90]);
            
            % Hack to avoid NaN issue because of shuffling
            dx = xi(2) - xi(1);
            edges = -90:dx:90;
            if edges(end) < 90
                edges = [edges, edges(end)+dx]; % extend so ≥90
            end
            
            % compute KDE again
            [f, xi] = ksdensity(grpOriErr_shuffled, edges, 'Function', 'pdf', 'Support', [-90 90]); % 'Support', [-90 90]
            
            % disp(trapz(xi, f))
            %f = f ./ trapz(xi, f);
            
            % Calculate likelihood on actual unshuffled data
            likelihood = calculateLikelihood(f, xi, grpOriErr);
            %assert(sum(isnan(likelihood)) == 0)
            % likelihood = likelihood.*pConfReport;        %
            likelihood = likelihood.*pConfReport_shuffled; % Model here is fit on shuffled data so pConf should come from shuffled data
            nll_ = -log(likelihood);
            totalNLL = totalNLL + sum(nll_);

        end
    end
end

fprintf('NLL by ori (shuffled conf label and uncertainty level): %.4f \n', totalNLL);


%% Distribution - Do only by ori
nitr = 2;
worstCaseNLLs = zeros(1 ,nitr);

for itr = 1:nitr
    
    confidence_report_all_shuffled = confidence_report_all( ...
        randperm(numel(confidence_report_all)));
    confidence_report_all_shuffled = reshape( ...
        confidence_report_all_shuffled, ...
        size(confidence_report_all));
    
    resp_err_all_shuffled = resp_err_all( ...
        randperm(numel(resp_err_all)));
    resp_err_all_shuffled = reshape( ...
        resp_err_all_shuffled, ...
        size(resp_err_all));
    
    % % Sanity test: equivalent to best case
    %confidence_report_all_shuffled = confidence_report_all;
    %resp_err_all_shuffled = resp_err_all;
    
    totalNLL = 0;
    
    for confVal = [0, 1]
        for i=1:n_uncertainty_levels
            for j=1:numel(orientations)
    
                oriErr_shuffled  = squeeze( resp_err_all_shuffled(i, j, :) );
                oriErr           = squeeze( resp_err_all(i, j, :) );
                reportedConf     = squeeze( confidence_report_all(i, j, :) );
                reportedConf_shuffled = squeeze( confidence_report_all_shuffled(i, j, :) );
                grpOriErr        = oriErr(reportedConf == confVal);
                %grpOriErr_shuffled = oriErr(reportedConf_shuffled == confVal);
                grpOriErr_shuffled = oriErr_shuffled(reportedConf_shuffled == confVal);
                
                pConfReport = sum(reportedConf == confVal) / numel(reportedConf);
                pConfReport_shuffled = sum(reportedConf_shuffled == confVal) / numel(reportedConf_shuffled);
                %pStim = 1/numel(orientations);
                
                % P(err, conf, stim) = P(err|conf,stim)P(conf|stim)P(stim)
                
                if numel(grpOriErr) == 0
                    continue
                end
                
                %assert(numel(grpOriErr) > 0);
                
                % Fit PDF on shuffled data
                [~, xi] = ksdensity(grpOriErr_shuffled, 'Function', 'pdf', 'Support', [-90 90]);
                
                % Hack to avoid NaN issue because of shuffling
                dx = xi(2) - xi(1);
                edges = -90:dx:90;
                if edges(end) < 90
                    edges = [edges, edges(end)+dx]; % extend so ≥90
                end
                
                % compute KDE again
                [f, xi] = ksdensity(grpOriErr_shuffled, edges, 'Function', 'pdf', 'Support', [-90 90]); % 'Support', [-90 90]
                
                % disp(trapz(xi, f))
                %f = f ./ trapz(xi, f);
                
                % Calculate likelihood on actual unshuffled data
                likelihood = calculateLikelihood(f, xi, grpOriErr);
                %assert(sum(isnan(likelihood)) == 0)
                % likelihood = likelihood.*pConfReport;        %
                likelihood = likelihood.*pConfReport_shuffled; % Model here is fit on shuffled data so pConf should come from shuffled data
                nll_ = -log(likelihood);
                totalNLL = totalNLL + sum(nll_);
    
            end
        end
    end
    
    fprintf('NLL by ori (shuffled conf label and uncertainty level): %.4f \n', totalNLL);

    worstCaseNLLs(itr) = totalNLL;

end

figure
histogram(worstCaseNLLs, DisplayName="Shuffled uncert and conf labels (FULL)")
hold on
xline(NLL_uniform_by_stim, LineStyle="--", DisplayName="All random")
xlabel("NLL")
ylabel("count")

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

%%
function likelihood = calculateLikelihood(pdf, x_pts, trlData_errors)
    likelihood = interp1(x_pts, pdf, trlData_errors, 'linear'); % eps 1e-10
end
