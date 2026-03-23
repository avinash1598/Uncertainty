function plotFitResult_guessrate(data, modelParams, initCond, modelType, errBins, fullModel)

if nargin < 6 && isempty(fullModel)
    fullModel = false;
end

grpOriErr            = data.resp_err_all; 
confReport           = data.confidence_report_all;
n_uncertainty_levels = size(grpOriErr, 1);
orientations         = data.orientations';

grpOriErr    = reshape(grpOriErr, n_uncertainty_levels, []);
confReport   = reshape(confReport, n_uncertainty_levels, []);
% errBins      = -90:3:90;
% errBins      = -90:0.1:90;

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

    if fullModel
        param_sigma_ori_scale = modelParams(n_uncertainty_levels + 5);
        % param_sigma_a    = modelParams(n_uncertainty_levels + 5:2*n_uncertainty_levels + 4);
        param_bias       = modelParams(n_uncertainty_levels + 6);
    end

elseif modelType == "ind"
    param_sigma_s        = modelParams(1:n_uncertainty_levels);
    param_shape          = modelParams(n_uncertainty_levels + 1);
    param_scale          = modelParams(n_uncertainty_levels + 2);
    param_sigma_meta     = modelParams(n_uncertainty_levels + 3);
    param_Cc             = modelParams(n_uncertainty_levels + 4);
    param_guessrate      = modelParams(n_uncertainty_levels + 5);

    if fullModel
        param_sigma_ori_scale = modelParams(n_uncertainty_levels + 6);
        % param_sigma_a    = modelParams(n_uncertainty_levels + 6:2*n_uncertainty_levels + 5);
        param_bias       = modelParams(n_uncertainty_levels + 7);
    end

else
    error("Invalid modelType")
end

analyticalSols = cell(n_uncertainty_levels, 1);
anlytcl_sigma_m = zeros(1, n_uncertainty_levels);
anlytcl_sigma_m_HC = zeros(1, n_uncertainty_levels);
anlytcl_sigma_m_LC = zeros(1, n_uncertainty_levels);
anlytcl_sigma_m_stim = zeros(n_uncertainty_levels, numel(orientations));
anlytcl_mad_m = zeros(1, n_uncertainty_levels);
anlytcl_mad_m_HC = zeros(1, n_uncertainty_levels);
anlytcl_mad_m_LC = zeros(1, n_uncertainty_levels);
anlytcl_mad_m_stim = zeros(n_uncertainty_levels, numel(orientations));
anlytcl_bias = zeros(n_uncertainty_levels, numel(orientations));

for i=1:n_uncertainty_levels

    mP.sigma_s             = param_sigma_s(i);
    mP.scale               = param_scale;
    mP.Cc                  = param_Cc;
    mP.sigma_meta          = param_sigma_meta;
    mP.guessRate           = param_guessrate;

    mP.sigma_m_shape1      = initCond.sigma_m_shape1;
    mP.sigma_m_shape2      = initCond.sigma_m_shape2;
    mP.biasShape           = initCond.biasShape;

    if fullModel
        mP.b                   = param_sigma_s(i);
        mP.a                   = param_sigma_ori_scale*param_sigma_s(i);
        % mP.a                   = param_sigma_a(i);
        mP.biasAmp             = param_bias;
    end
    
    if modelType == "ind"
        mP.shape = param_shape;

        if fullModel
            analyticalSol = getEstimatesPDFs(orientations, mP, false);
        else
            analyticalSol = getEstimatesPDFs_reduced_model(errBins, mP, false);
        end
    
    elseif modelType == "cov"
        
        if fullModel
            analyticalSol = getEstimationsPDF_cov(orientations, mP, false);
        else
            analyticalSol = getEstimationsPDF_cov_reduced(errBins, mP, false);
        end
    end
    
    anlytcl_sigma_m(i)    = analyticalSol.E_sigma_m;
    anlytcl_sigma_m_HC(i) = analyticalSol.E_sigma_m_HC;
    anlytcl_sigma_m_LC(i) = analyticalSol.E_sigma_m_LC;

    anlytcl_mad_m(i)         = analyticalSol.mad_m;
    anlytcl_mad_m_HC(i)      = analyticalSol.mad_m_HC;
    anlytcl_mad_m_LC(i)      = analyticalSol.mad_m_LC;
    
    if fullModel
        anlytcl_mad_m_stim(i, :) = analyticalSol.mad_m_by_ori;
        anlytcl_bias(i, :)       = analyticalSol.bias;
        anlytcl_sigma_m_stim(i, :) = analyticalSol.E_sigma_m_stim;
    end
    
    analyticalSols{i}     = analyticalSol;
end

% Is this right??
analytical_bias = analyticalSol.bias;

warning("only for full model. Change it for reduced as well")
% Arrange in increasing order of sigma_s
sigma_s_stim = param_sigma_s' + ...
    param_sigma_ori_scale.*param_sigma_s'*( ...
    abs(sind(initCond.sigma_m_shape1*orientations - initCond.sigma_m_shape2)));
sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(analytical_bias).^2 ); % wrong % param_bias
% sigma_s_reduced_model = sqrt( mean( sigma_s_stim.^2, 2 ) + std(param_bias).^2 )';

[~, idxSorted] = sort(sigma_s_reduced_model);

% Plot PDFs
figure 

for i=1:n_uncertainty_levels
    
    idx = idxSorted(i);
    
    cR = confReport(idx, :);
    dataHC = grpOriErr(idx, cR == 1);
    dataLC = grpOriErr(idx, cR == 0);
    
    subplot(2, n_uncertainty_levels, i);
    hold on;
    
    % y = pdf_stim_LC(idx, :);
    % scatter(errBins, y(:), 'filled', DisplayName="LC");
    histogram(dataLC, errBins, Normalization="pdf");
    plot(analyticalSols{idx}.rvOriErrs, analyticalSols{idx}.analyticalPDF_LC, DisplayName="fit", LineWidth=1.5);
    hold off
    
    %xline(0, LineStyle="--"); Probalamatic line
    xlabel("Error (deg)");
    ylabel("P( Err / LC )");
    title(sprintf("C: %.2f, S: %d, D: %.2f",  ...
        data.uncertaintyVals(idx, 1), ...
        data.uncertaintyVals(idx, 2), ...
        data.uncertaintyVals(idx, 3)))
    
    legend;
    
    subplot(2, n_uncertainty_levels, n_uncertainty_levels + i);
    hold on;

    % y = pdf_stim_HC(idx, :);
    % scatter(errBins, y(:), 'filled', DisplayName="HC");
    histogram(dataHC, errBins, Normalization="pdf");
    plot(analyticalSols{idx}.rvOriErrs, analyticalSols{idx}.analyticalPDF_HC, DisplayName="fit", LineWidth=1.5);
    
    hold off
    
    %xline(0, LineStyle="--"); Probalamatic line
    xlabel("Error (deg)");
    ylabel("P( Err / HC )");

    legend;
    hold off;

end

figure 

for i=1:n_uncertainty_levels

    idx = idxSorted(i);
    errData = grpOriErr(idx, :);

    subplot(2, n_uncertainty_levels/2, i);
    hold on;

    % y = pdf_stim_HC(idx, :);
    % scatter(errBins, y(:), 'filled', DisplayName="HC");
    histogram(errData, errBins, Normalization="pdf");
    plot(analyticalSols{idx}.rvOriErrs, analyticalSols{idx}.analyticalPDF, DisplayName="fit", LineWidth=1.5);
    
    hold off
    
    %xline(0, LineStyle="--"); Probalamatic line
    xlabel("Error (deg)");
    ylabel("P( Err )");
    title(sprintf("C: %.2f, S: %d, D: %.2f",  ...
        data.uncertaintyVals(idx, 1), ...
        data.uncertaintyVals(idx, 2), ...
        data.uncertaintyVals(idx, 3)))

    legend;
    hold off;
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
errorbar(sigma_s_reduced_model(idxSorted), ...
    x(idxSorted), y(idxSorted), ...
    'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName="High confidence");

xlabel("\sigma_s")
% xticks(round( param_sigma_s(idxSorted), 1 ))
% xticklabels(round( param_sigma_s(idxSorted), 1 ))
ylabel("Error")

subplot(2, 3, 2)

% Behavioral variability
scatter(sigma_s_reduced_model(idxSorted), y(idxSorted), "filled");
hold on
plot(sigma_s_reduced_model(idxSorted), anlytcl_sigma_m(idxSorted), LineWidth=1.5);
xlabel("\sigma_s (sensory noise)")
ylabel("\sigma_m (measurement noise)")
% xticks(round( param_sigma_s(idxSorted), 1 ))
% xticklabels(round( param_sigma_s(idxSorted), 1 ))
hold off

subplot(2, 3, 3)


% Behavioral variability
scatter(sigma_s_reduced_model(idxSorted), y_HC(idxSorted), "filled", DisplayName="High confidence");
hold on
plot(sigma_s_reduced_model(idxSorted), anlytcl_sigma_m_HC(idxSorted), LineWidth=1.5, HandleVisibility="off");
scatter(sigma_s_reduced_model(idxSorted), y_LC(idxSorted), "filled", DisplayName="Low confidence");
plot(sigma_s_reduced_model(idxSorted), anlytcl_sigma_m_LC(idxSorted), LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s (sensory noise)")
ylabel("\sigma_m (measurement noise)")
% xticks(round( param_sigma_s(idxSorted), 1 ))
% xticklabels(round( param_sigma_s(idxSorted), 1 ))
legend
hold off

subplot(2, 3, 4)

% Behavioral variability
scatter(sigma_s_reduced_model(idxSorted).^2, y_HC(idxSorted).^2, "filled", DisplayName="High confidence");
hold on
plot(sigma_s_reduced_model(idxSorted).^2, anlytcl_sigma_m_HC(idxSorted).^2, LineWidth=1.5, HandleVisibility="off");
scatter(sigma_s_reduced_model(idxSorted).^2, y_LC(idxSorted).^2, "filled", DisplayName="Low confidence");
plot(sigma_s_reduced_model(idxSorted).^2, anlytcl_sigma_m_LC(idxSorted).^2, LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s^2 (sensory noise)")
ylabel("\sigma_m^2 (measurement noise)")
% xticks(round( param_sigma_s(idxSorted), 1 ).^2)
% xticklabels(round( param_sigma_s(idxSorted), 1 ).^2)
legend
hold off

subplot(2, 3, 5)

% Behavioral variability
scatter(sigma_s_reduced_model(idxSorted), y_m(idxSorted), "filled");
hold on
plot(sigma_s_reduced_model(idxSorted), anlytcl_mad_m(idxSorted), LineWidth=1.5);
xlabel("\sigma_s (sensory noise)")
ylabel("MAD (measurement)")
hold off

subplot(2, 3, 6)

% Behavioral variability
scatter(sigma_s_reduced_model(idxSorted), y_HC_m(idxSorted), "filled", DisplayName="High confidence");
hold on
plot(sigma_s_reduced_model(idxSorted), anlytcl_mad_m_HC(idxSorted), LineWidth=1.5, HandleVisibility="off");
scatter(sigma_s_reduced_model(idxSorted), y_LC_m(idxSorted), "filled", DisplayName="Low confidence");
plot(sigma_s_reduced_model(idxSorted), anlytcl_mad_m_LC(idxSorted), LineWidth=1.5, HandleVisibility="off");
xlabel("\sigma_s(s) (sensory noise)")
ylabel("MAD (measurement)")
legend
hold off

if fullModel

    % Ori dependent variance
    figure

    stdByOri = data.stdByOri;
    % madByOri = data.madByOri;
    % bias = data.bias;
    orientations = data.orientations;
    %b_opt = modelParams(1:n_uncertainty_levels);
    %a_opt = modelParams(end-1).*b_opt;
    
    for i=1:n_uncertainty_levels
    
        %b_est = b_opt(i);
        %a_est = a_opt(i);
        
        subplot(2, 3, i)
        hold on
        %plot( orientations, b_est + a_est*abs(sind(orientations - 90)), LineWidth=1.5, DisplayName="fit")
        %plot( orientations, b_est + a_est*abs(sind(2*orientations)), LineWidth=1.5, DisplayName="fit")
        scatter(orientations, stdByOri(i, :))
        plot(orientations, anlytcl_sigma_m_stim(i, :), LineWidth=1.5, DisplayName="fit")
        % scatter(orientations, madByOri(i, :))
        % ylim([0 60])
        hold off
    end
    
    % Bias
    figure
    
    bias = data.bias;
    oriBias = anlytcl_bias(1, :);
    %biasAmp = modelParams(end);
    %oriBias = biasAmp*sind(2*orientations);
    
    scatter(orientations, bias)
    hold on
    plot(orientations, oriBias)
    hold off

end

end
