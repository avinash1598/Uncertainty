function retData = formatExpData(data, fltData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRT actual stim orientation
% TODO: Also return bias corrected bhv errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimOris       = data.dat.stimOri;
orientations   = unique(stimOris);
n_orientations = numel(orientations);

% Error wrt sample mean
sampleMeanOri  = data.dat.stimSampleMeanOri;
reportedOri    = data.dat.reportedOri;

rawOriError_S = reportedOri - sampleMeanOri;
rawOriError_S = mod(rawOriError_S + 90, 180) - 90;

% Error wrt actual orientation
rawOriError = data.dat.rawOriError;

[G, c, s, d]         = findgroups(data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur);
grpIdxes             = unique(G);
n_uncertainty_levels = numel(grpIdxes);
nTrials              = numel(rawOriError(G == grpIdxes(1) & (stimOris == orientations(1))));

uncertaintyVals = [c s d];

% Arrange groups in increasing uncertainty level
fitStds = zeros(1, n_uncertainty_levels);

for i=1:n_uncertainty_levels
    grpIdx          = grpIdxes(i);
    grpRawErr_Flt   = rawOriError(G == grpIdx);
    
    pd = fitdist(grpRawErr_Flt(~isnan(grpRawErr_Flt)), 'Normal'); % there should be non nan so need to filter nans
    
    mu = pd.mu;
    sigma = pd.sigma;
    fitStds(i) = sigma;
end

[~, idx] = sort(fitStds);
sidx = grpIdxes(idx);
uncertaintyVals = uncertaintyVals(sidx, :);

% Arrange in structured format - same as the one used for analysis
theta_true_all        = zeros(n_uncertainty_levels, n_orientations, nTrials);
theta_resp_all        = zeros(n_uncertainty_levels, n_orientations, nTrials); % Recorded theta based on user response
confidence_report_all = zeros(n_uncertainty_levels, n_orientations, nTrials);
resp_err_all          = zeros(n_uncertainty_levels, n_orientations, nTrials);
theta_true_all_S      = zeros(n_uncertainty_levels, n_orientations, nTrials);
resp_err_all_S        = zeros(n_uncertainty_levels, n_orientations, nTrials);

stim_contrast_all     = zeros(n_uncertainty_levels, n_orientations, nTrials);
stim_dispersion_all   = zeros(n_uncertainty_levels, n_orientations, nTrials);
stim_dur_all          = zeros(n_uncertainty_levels, n_orientations, nTrials);

for i=1:n_uncertainty_levels
    for j = 1:n_orientations
        
        grpIdx = sidx(i);
        fltIDx = ( G == grpIdx) & (stimOris == orientations(j)) ;

        % WRT actual stim ori
        grpRawErr_Flt   = rawOriError(fltIDx);
        grpStimOri      = data.dat.stimOri(fltIDx);
        grpOriResp      = data.dat.reportedOri(fltIDx);
        grpReportedConf = data.dat.reportedConf(fltIDx);
        
        theta_true_all(i, j, :)          = grpStimOri(:);
        theta_resp_all(i, j, :)          = grpOriResp(:);
        confidence_report_all(i, j, :)   = grpReportedConf(:);
        resp_err_all(i, j, :)            = grpRawErr_Flt(:);
        
        stim_contrast_all(i, j, :)       = data.dat.stimContrast(fltIDx);
        stim_dispersion_all(i, j, :)     = data.dat.stimSpread(fltIDx);
        stim_dur_all(i, j, :)            = data.dat.stimDur(fltIDx);
        
        % WRT sample mean
        grpRawErr_Flt_S   = rawOriError_S(fltIDx);
        grpStimOri_S      = data.dat.stimSampleMeanOri(fltIDx);
        
        theta_true_all_S(i, j, :)          = grpStimOri_S(:);
        resp_err_all_S(i, j, :)            = grpRawErr_Flt_S(:);

    end
end

retData.uncertaintyVals       = uncertaintyVals;
retData.n_uncertainty_levels  = n_uncertainty_levels;
retData.theta_true_all        = theta_true_all;
retData.theta_resp_all        = theta_resp_all;
retData.confidence_report_all = confidence_report_all;
retData.resp_err_all          = resp_err_all;
retData.theta_true_all_S      = theta_true_all_S;
retData.resp_err_all_S        = resp_err_all_S;
retData.resp_err_all_S        = resp_err_all_S;

retData.stim_contrast_all     = stim_contrast_all;
retData.stim_dispersion_all   = stim_dispersion_all;
retData.stim_dur_all          = stim_dur_all;

end
