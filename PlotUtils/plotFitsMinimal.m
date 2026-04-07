function plotFitsMinimal(data, modelParams, initCond, modelType, fitType, pltIdx, plotMADs)

assert( fitType == "full" );

if nargin < 7
    plotMADs = false;
end

grpOriErr            = data.resp_err_all; 
confReport           = data.confidence_report_all;
n_uncertainty_levels = size(grpOriErr, 1);
orientations         = data.orientations';

grpOriErr    = reshape(grpOriErr, n_uncertainty_levels, []);
confReport   = reshape(confReport, n_uncertainty_levels, []);

% Analytical solution
if modelType == "cov"
    param_sigma_s        = modelParams(1:n_uncertainty_levels);
    param_scale          = modelParams(n_uncertainty_levels + 1);
    param_sigma_meta     = modelParams(n_uncertainty_levels + 2);
    param_Cc             = modelParams(n_uncertainty_levels + 3);
    param_exp            = modelParams(n_uncertainty_levels + 4);
    param_guessrate      = modelParams(n_uncertainty_levels + 5);
    
    if fitType == "full"
        param_sigma_ori_scale = modelParams(n_uncertainty_levels + 6);
        param_bias       = modelParams(n_uncertainty_levels + 7);
    end

elseif modelType == "ind"
    param_sigma_s        = modelParams(1:n_uncertainty_levels);
    param_shape          = modelParams(n_uncertainty_levels + 1);
    param_scale          = modelParams(n_uncertainty_levels + 2);
    param_sigma_meta     = modelParams(n_uncertainty_levels + 3);
    param_Cc             = modelParams(n_uncertainty_levels + 4);
    param_guessrate      = modelParams(n_uncertainty_levels + 5);

    if fitType == "full"
        param_sigma_ori_scale = modelParams(n_uncertainty_levels + 6);
        param_bias       = modelParams(n_uncertainty_levels + 7);
    end

elseif modelType == "singlyStochastic"
    param_sigma_s        = modelParams(1:n_uncertainty_levels);
    param_sigma_meta     = modelParams(n_uncertainty_levels + 1);
    param_Cc             = modelParams(n_uncertainty_levels + 2);
    param_guessrate      = modelParams(n_uncertainty_levels + 3);
    
    if fitType == "full"
        param_sigma_ori_scale = modelParams(n_uncertainty_levels + 4);
        param_bias            = modelParams(n_uncertainty_levels + 5);
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
anlytcl_sigma_s_reduced = zeros(1, n_uncertainty_levels);

for i=1:n_uncertainty_levels

    mP.sigma_s             = param_sigma_s(i);
    mP.Cc                  = param_Cc;
    mP.sigma_meta          = param_sigma_meta;
    mP.guessRate           = param_guessrate;

    mP.sigma_m_shape1      = initCond.sigma_m_shape1;
    mP.sigma_m_shape2      = initCond.sigma_m_shape2;
    mP.biasShape           = initCond.biasShape;

    if fitType == "full"
        mP.b                   = param_sigma_s(i);
        mP.a                   = param_sigma_ori_scale*param_sigma_s(i);
        mP.biasAmp             = param_bias;
    end
    
    if modelType == "ind"
        mP.scale = param_scale;
        mP.shape = param_shape;
        
        if fitType == "full"
            analyticalSol = getPDFs_ind(orientations, mP, false);
        else
            analyticalSol = getPDFs_ind_reduced(mP, false);
        end
    
    elseif modelType == "cov"
        mP.scale = param_scale;
        mP.exp   = param_exp;

        if fitType == "full"
            analyticalSol = getPDFs_cov(orientations, mP, false);
        else
            analyticalSol = getPDFs_cov_reduced(mP, false);
        end

    elseif modelType == "singlyStochastic"
        if fitType == "full"
            analyticalSol = getPDFs_SinglyStochastic(orientations, mP, false);
        else
            error("No reduced version of this model!")
        end
    end
    
    anlytcl_sigma_m(i)    = analyticalSol.E_sigma_m;
    anlytcl_sigma_m_HC(i) = analyticalSol.E_sigma_m_HC;
    anlytcl_sigma_m_LC(i) = analyticalSol.E_sigma_m_LC;

    anlytcl_mad_m(i)         = analyticalSol.mad_m;
    anlytcl_mad_m_HC(i)      = analyticalSol.mad_m_HC;
    anlytcl_mad_m_LC(i)      = analyticalSol.mad_m_LC;

    anlytcl_sigma_s_reduced(i) = analyticalSol.analytical_sigma_s_reduced;

    if fitType == "full" || fitType == "full_jumbo"
        anlytcl_mad_m_stim(i, :) = analyticalSol.mad_m_by_ori;
        anlytcl_bias(i, :)       = analyticalSol.bias;
        anlytcl_sigma_m_stim(i, :) = analyticalSol.E_sigma_m_stim;
    end
    
    analyticalSols{i}     = analyticalSol;
end

% Is this right??
sigma_s_reduced_model = anlytcl_sigma_s_reduced;

[~, idxSorted] = sort(sigma_s_reduced_model);

%% Plot results

x = mean(grpOriErr, 2);
y = std(grpOriErr, 0, 2);
y_m = mad(grpOriErr, 1, 2);

HC_idx = confReport == 1;
LC_idx = confReport == 0;

resp_HC = grpOriErr;
resp_HC(~HC_idx) = NaN;

resp_LC = grpOriErr;
resp_LC(~LC_idx) = NaN;

y_HC = std(resp_HC, 0, 2, 'omitnan');
y_HC_m = mad(resp_HC, 1, 2);

y_LC = std(resp_LC, 0, 2, 'omitnan');
y_LC_m = mad(resp_LC, 1, 2);

if ~plotMADs

    subplot(pltIdx(1), pltIdx(2), pltIdx(3))
    
    % Behavioral variability
    scatter(sigma_s_reduced_model(idxSorted), y(idxSorted), "filled");
    hold on
    plot(sigma_s_reduced_model(idxSorted), anlytcl_sigma_m(idxSorted), LineWidth=1.5);
    xlabel("\sigma_s (sensory noise)")
    ylabel("\sigma_m (measurement noise)")
    hold off

    subplot(pltIdx(1), pltIdx(2), pltIdx(3)+1)
    
    % Behavioral variability
    scatter(sigma_s_reduced_model(idxSorted), y_HC(idxSorted), "filled", DisplayName="High confidence");
    hold on
    plot(sigma_s_reduced_model(idxSorted), anlytcl_sigma_m_HC(idxSorted), LineWidth=1.5, HandleVisibility="off");
    scatter(sigma_s_reduced_model(idxSorted), y_LC(idxSorted), "filled", DisplayName="Low confidence");
    plot(sigma_s_reduced_model(idxSorted), anlytcl_sigma_m_LC(idxSorted), LineWidth=1.5, HandleVisibility="off");
    xlabel("\sigma_s (sensory noise)")
    ylabel("\sigma_m (measurement noise)")
    legend
    hold off

else
    subplot(pltIdx(1), pltIdx(2), pltIdx(3))

    % Behavioral variability
    scatter(sigma_s_reduced_model(idxSorted), y_m(idxSorted), "filled");
    hold on
    plot(sigma_s_reduced_model(idxSorted), anlytcl_mad_m(idxSorted), LineWidth=1.5);
    xlabel("\sigma_s (sensory noise)")
    ylabel("MAD (measurement)")
    hold off
    
    subplot(pltIdx(1), pltIdx(2), pltIdx(3)+1)
    
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

end

end
