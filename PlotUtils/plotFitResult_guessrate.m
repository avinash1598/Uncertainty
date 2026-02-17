function plotFitResult_guessrate(data, modelParams, modelType)

addpath('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Scripts/CompleteModelEstimation/LL_scripts/')

grpOriErr            = data.resp_err_all; 
confReport           = data.confidence_report_all;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr    = reshape(grpOriErr, n_uncertainty_levels, []);
confReport   = reshape(confReport, n_uncertainty_levels, []);
errBins      = -90:3:90;

% Get PDFs from data for HC and LC
pdf_stim_LC = zeros( n_uncertainty_levels, numel(errBins) );
pdf_stim_HC = zeros( n_uncertainty_levels, numel(errBins) );

for i=1:n_uncertainty_levels

    cR = confReport(i, :);
    dataHC = grpOriErr(i, cR == 1);
    dataLC = grpOriErr(i, cR == 0);
    
    centers = errBins;
    binWidth = mean(diff(centers));
    edges = [centers - binWidth/2, centers(end) + binWidth/2];
    
    [pdfHC, ~] = histcounts(dataHC, ...
        'Normalization', 'pdf', ...
        'BinEdges', edges);

    [pdfLC, ~] = histcounts(dataLC, ...
        'Normalization', 'pdf', ...
        'BinEdges', edges);
    
    pdf_stim_HC(i, :) = pdfHC;
    pdf_stim_LC(i, :) = pdfLC;
    
end

% Analytical solution
if modelType == "cov"
    param_sigma_s        = modelParams(1:n_uncertainty_levels);
    param_scale          = modelParams(n_uncertainty_levels + 1);
    param_sigma_meta     = modelParams(n_uncertainty_levels + 2);
    param_Cc             = modelParams(n_uncertainty_levels + 3);
    param_guessrate      = modelParams(n_uncertainty_levels + 4);

elseif modelType == "ind"
    param_sigma_s        = modelParams(1:n_uncertainty_levels);
    param_shape          = modelParams(n_uncertainty_levels + 1);
    param_scale          = modelParams(n_uncertainty_levels + 2);
    param_sigma_meta     = modelParams(n_uncertainty_levels + 3);
    param_Cc             = modelParams(n_uncertainty_levels + 4);
    param_guessrate      = modelParams(n_uncertainty_levels + 5);
else
    error("Invalid modelType")
end

analyticalSols = cell(n_uncertainty_levels, 1);
anlytcl_sigma_m = zeros(1, n_uncertainty_levels);
anlytcl_sigma_m_HC = zeros(1, n_uncertainty_levels);
anlytcl_sigma_m_LC = zeros(1, n_uncertainty_levels);
anlytcl_mad_m_stim = zeros(1, n_uncertainty_levels);
anlytcl_mad_m_stim_HC = zeros(1, n_uncertainty_levels);
anlytcl_mad_m_stim_LC = zeros(1, n_uncertainty_levels);


for i=1:n_uncertainty_levels

    mP.sigma_s             = param_sigma_s(i);
    mP.scale               = param_scale;
    mP.Cc                  = param_Cc;
    mP.sigma_meta          = param_sigma_meta;
    mP.guessRate           = param_guessrate;

    if modelType == "ind"
        mP.shape = param_shape;
        analyticalSol = getEstimatesPDFs_reduced_model(errBins, mP);
    elseif modelType == "cov"
        analyticalSol = getEstimationsPDF_cov_reduced(errBins, mP);
    end
    
    anlytcl_sigma_m(i)    = analyticalSol.E_sigma_m;
    anlytcl_sigma_m_HC(i) = analyticalSol.E_sigma_m_HC;
    anlytcl_sigma_m_LC(i) = analyticalSol.E_sigma_m_LC;

    anlytcl_mad_m_stim(i)      = analyticalSol.mad_m;
    anlytcl_mad_m_stim_HC(i)   = analyticalSol.mad_m_HC;
    anlytcl_mad_m_stim_LC(i)   = analyticalSol.mad_m_LC;

    analyticalSols{i}     = analyticalSol;
end


% Arrange in increasing order of sigma_s
[~, idxSorted] = sort(param_sigma_s);

% Plot PDFs
figure 

for i=1:n_uncertainty_levels

    idx = idxSorted(i);

    subplot(2, n_uncertainty_levels/2, i)
    hold on
    
    y = pdf_stim_LC(idx, :);
    scatter(errBins, y(:), 'filled', DisplayName="LC");
    plot(errBins, analyticalSols{idx}.analyticalPDF_LC, HandleVisibility="off", LineWidth=1.5)
    
    y = pdf_stim_HC(idx, :);
    scatter(errBins, y(:), 'filled', DisplayName="HC");
    plot(errBins, analyticalSols{idx}.analyticalPDF_HC, HandleVisibility="off", LineWidth=1.5)
    
    xline(0, LineStyle="--")
    xlabel("Error (deg)")
    ylabel("count")
    % title("")

    legend
    hold off

end

% Plot error
%% Plot results

figure

x = mean(grpOriErr, 2);
y = std(grpOriErr, 0, 2);
y_m = mad(grpOriErr, 1, 2);

HC_idx = confReport == 1;
LC_idx = confReport == 0;

resp_HC = grpOriErr;
resp_HC(~HC_idx) = NaN;

resp_LC = grpOriErr;
resp_LC(~LC_idx) = NaN;

% x_HC = mean(resp_HC, 2, 'omitnan');
y_HC = std(resp_HC, 0, 2, 'omitnan');
y_HC_m = mad(resp_HC, 1, 2);

% x_LC = mean(resp_LC, 2, 'omitnan');
y_LC = std(resp_LC, 0, 2, 'omitnan');
y_LC_m = mad(resp_LC, 1, 2);

% x1 = resp_HC(1, :); valid_idx = ~isnan(x1); x1 = x1(valid_idx);
% x2 = resp_LC(1, :); valid_idx = ~isnan(x2); x2 = x2(valid_idx);
% x3 = resp_HC(n_uncertainty_levels, :); valid_idx = ~isnan(x3); x3 = x3(valid_idx);
% x4 = resp_LC(n_uncertainty_levels, :); valid_idx = ~isnan(x4); x4 = x4(valid_idx);

subplot(2, 3, 1)
errorbar(param_sigma_s(idxSorted), ...
    x(idxSorted), y(idxSorted), ...
    'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("\sigma_s")
% xticks(round( param_sigma_s(idxSorted), 1 ))
% xticklabels(round( param_sigma_s(idxSorted), 1 ))
ylabel("Error")

subplot(2, 3, 2)

% Behavioral variability
scatter(param_sigma_s(idxSorted), y(idxSorted), "filled");
hold on
plot(param_sigma_s(idxSorted), anlytcl_sigma_m(idxSorted), LineWidth=1.5);
xlabel("\sigma_s (sensory noise)")
ylabel("\sigma_m (measurement noise)")
% xticks(round( param_sigma_s(idxSorted), 1 ))
% xticklabels(round( param_sigma_s(idxSorted), 1 ))
hold off

subplot(2, 3, 3)


% Behavioral variability
scatter(param_sigma_s(idxSorted), y_HC(idxSorted), "filled", DisplayName="High confidence");
hold on
plot(param_sigma_s(idxSorted), anlytcl_sigma_m_HC(idxSorted), LineWidth=1.5, HandleVisibility="off");
scatter(param_sigma_s(idxSorted), y_LC(idxSorted), "filled", DisplayName="Low confidence");
plot(param_sigma_s(idxSorted), anlytcl_sigma_m_LC(idxSorted), LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s (sensory noise)")
ylabel("\sigma_m (measurement noise)")
% xticks(round( param_sigma_s(idxSorted), 1 ))
% xticklabels(round( param_sigma_s(idxSorted), 1 ))
legend
hold off

subplot(2, 3, 4)

% Behavioral variability
scatter(param_sigma_s(idxSorted).^2, y_HC(idxSorted).^2, "filled", DisplayName="High confidence");
hold on
plot(param_sigma_s(idxSorted).^2, anlytcl_sigma_m_HC(idxSorted).^2, LineWidth=1.5, HandleVisibility="off");
scatter(param_sigma_s(idxSorted).^2, y_LC(idxSorted).^2, "filled", DisplayName="Low confidence");
plot(param_sigma_s(idxSorted).^2, anlytcl_sigma_m_LC(idxSorted).^2, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s^2 (sensory noise)")
ylabel("\sigma_m^2 (measurement noise)")
% xticks(round( param_sigma_s(idxSorted), 1 ).^2)
% xticklabels(round( param_sigma_s(idxSorted), 1 ).^2)
legend
hold off

subplot(2, 3, 5)

% Behavioral variability
scatter(param_sigma_s(idxSorted), y_m(idxSorted), "filled");
hold on
plot(param_sigma_s(idxSorted), anlytcl_mad_m_stim(idxSorted), LineWidth=1.5);
xlabel("\sigma_s (sensory noise)")
ylabel("MAD (measurement)")
hold off

subplot(2, 3, 6)

% Behavioral variability
scatter(param_sigma_s(idxSorted), y_HC_m(idxSorted), "filled", DisplayName="High confidence");
hold on
plot(param_sigma_s(idxSorted), anlytcl_mad_m_stim_HC(idxSorted), LineWidth=1.5, HandleVisibility="off");
scatter(param_sigma_s(idxSorted), y_LC_m(idxSorted), "filled", DisplayName="Low confidence");
plot(param_sigma_s(idxSorted), anlytcl_mad_m_stim_LC(idxSorted), LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (sensory noise)")
ylabel("MAD (measurement)")
legend
hold off


end
