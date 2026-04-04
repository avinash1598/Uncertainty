function nlls = getNLLForEachTrial(data, initCond, modelType, fitType, modelParams, fltTrlIdx)
    assert(fitType == "reduced" || fitType == "full" || fitType == "full_jumbo");

    trlData              = convertToTrialData(data);
    % grpOriErr            = trlData.grpOriErr;
    n_uncertainty_levels = trlData.n_uncertainty_levels;
    trlErrors            = trlData.trlErrors;
    trlConfReports       = trlData.trlConfReports;
    trlUncertaintyLevels = trlData.trlUncertaintyLevels;
    trlStimOris          = trlData.trlStimOris;
     
    if nargin > 5 && ~isempty(fltTrlIdx)
        trlErrors            = trlErrors(fltTrlIdx);
        trlConfReports       = trlConfReports(fltTrlIdx);
        trlUncertaintyLevels = trlUncertaintyLevels(fltTrlIdx);
        trlStimOris          = trlStimOris(fltTrlIdx);
    end
        
    metaData.n_levels             = n_uncertainty_levels;
    metaData.orientations         = data.orientations';
    metaData.trlErrors            = trlErrors;
    metaData.trlConfReports       = trlConfReports;
    metaData.trlUncertaintyLevels = trlUncertaintyLevels;
    metaData.trlStimOris          = trlStimOris;
    metaData.sigma_m_shape1       = initCond.sigma_m_shape1; % 1
    metaData.sigma_m_shape2       = initCond.sigma_m_shape2; % 90
    metaData.biasShape            = initCond.biasShape;      % 2
    metaData.a                    = initCond.a;
    metaData.biasAmp              = initCond.biasAmp;
    
    if modelType == "cov"
        nlls = computeNLLCov(modelParams, metaData, fitType, false); 
    elseif modelType == "ind"
        nlls = computeNLL(modelParams, metaData, fitType, false); 
    else
        nlls = nan;
    end
end
