function retData = formatExpData(data, debias, sortByMAD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRT actual stim orientation
% TODO: Also return bias corrected bhv errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: remove error bias. Assuming bias does not depend upon uncertainty

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

if debias
    % debias
    stimSummary = groupsummary(data.dat, {'stimOri'}, {'mean', 'std'}, 'rawOriError');
    bias     = stimSummary.mean_rawOriError;
    stimVals = stimSummary.stimOri;
    
    % Error wrt to stim ori
    for i = 1:numel(stimVals)
        idx = stimOris == stimVals(i);
        rawOriError(idx) = rawOriError(idx) - bias(i);
    end
    
    % Error wrt sample mean
    for i = 1:numel(stimVals)
        idx = stimOris == stimVals(i);
        rawOriError_S(idx) = rawOriError_S(idx) - bias(i);
    end

end

[G, c, s, d]         = findgroups(data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur);
grpIdxes             = unique(G);
n_uncertainty_levels = numel(grpIdxes);
nTrials              = numel(rawOriError(G == grpIdxes(1) & (stimOris == orientations(1))));

uncertaintyVals = [c s d];

% Arrange groups in increasing uncertainty level (this might not be fully
% accurate or is it?) Maybe increasing sensory noise is right approach
fitStds = zeros(1, n_uncertainty_levels);
mads    = zeros(1, n_uncertainty_levels);

for i=1:n_uncertainty_levels
    grpIdx          = grpIdxes(i);
    grpRawErr_Flt   = rawOriError(G == grpIdx);
    
    assert(all(~isnan(grpRawErr_Flt))) % there should be no nan
    pd = fitdist(grpRawErr_Flt, 'Normal'); % there should be non nan so need to filter nans

    mu = pd.mu;
    sigma = pd.sigma;
    fitStds(i) = sigma;
    mads(i) = mad(grpRawErr_Flt, 1);
end

if sortByMAD
    [~, idx] = sort(mads);
else
    [~, idx] = sort(fitStds);
end

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


%% Marginal distributions
x = uncertaintyVals; % Size 3: contrast, spread, duration

assert(size(x, 2) == 3); % size != 3 => problem

combinations = {};

% i == 1: spread and duration
% i == 2: contrast and duration
% i == 3: contrast and spread
for i = 1:3 
    vals = x(:, i);

    [uniqVals, ~, idx] = unique(vals);  
    counts = accumarray(idx, 1);         
    
    indsUnique = find(counts == 4); % Find which unique values occur 4 times

    if numel(indsUnique) > 0
        valuesOccur4 = uniqVals(indsUnique); % Get the actual values that occur 4 times
        
        % Get all indices in vals where those values occur
        indicesInVals = vals == valuesOccur4; %find(ismember(vals, valuesOccur4));

        if i == 1
            combinations.sp_dur = x(indicesInVals, :, :);
        elseif i == 2
            combinations.ct_dur = x(indicesInVals, :, :);
        else
            combinations.ct_sp  = x(indicesInVals, :, :);
        end
    end

end

% Contrast and spread
tableVals      = [data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur]; % Combine relevant columns from the table
rowsToKeep     = ismember(tableVals, combinations.ct_sp, 'rows');                % Find rows matching any combination in factorialCombinations1
filteredTable  = data.dat(rowsToKeep, :); % Filter table

% 1: contrast groups ==============================================
[G, groupValues]     = findgroups(filteredTable.stimContrast);
grpIdxes             = unique(G);
n_uncertainty_levels = numel(grpIdxes);

sidx = getSortedGrpIdxes(n_uncertainty_levels, G, grpIdxes, filteredTable);

marginal_contrast_theta_true_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both contrast level no of trials must be equal
marginal_contrast_theta_resp_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both contrast level no of trials must be equal
marginal_contrast_resp_err_all          = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both contrast level no of trials must be equal
marginal_contrast_confidence_report_all = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both contrast level no of trials must be equal

marginal_contrast_vals = groupValues(sidx);

for i=1:n_uncertainty_levels
        
    grpIdx = sidx(i);
    fltIDx = ( G == grpIdx) ;
    
    % WRT actual stim ori
    grpRawErr_Flt   = filteredTable.rawOriError(fltIDx);
    grpStimOri      = filteredTable.stimOri(fltIDx);
    grpOriResp      = filteredTable.reportedOri(fltIDx);
    grpReportedConf = filteredTable.reportedConf(fltIDx);
    
    marginal_contrast_theta_true_all(i, :)          = grpStimOri(:);
    marginal_contrast_theta_resp_all(i, :)          = grpOriResp(:);
    marginal_contrast_confidence_report_all(i, :)   = grpReportedConf(:);
    marginal_contrast_resp_err_all(i, :)            = grpRawErr_Flt(:);
end

retData.marginal_contrast_vals                   = marginal_contrast_vals;
retData.marginal_contrast_theta_true_all         = marginal_contrast_theta_true_all;
retData.marginal_contrast_theta_resp_all         = marginal_contrast_theta_resp_all;
retData.marginal_contrast_resp_err_all           = marginal_contrast_resp_err_all;
retData.marginal_contrast_confidence_report_all  = marginal_contrast_confidence_report_all;


% 1: spread groups ==============================================
[G, groupValues]     = findgroups(filteredTable.stimSpread);
grpIdxes             = unique(G);
n_uncertainty_levels = numel(grpIdxes);

sidx = getSortedGrpIdxes(n_uncertainty_levels, G, grpIdxes, filteredTable);

marginal_spread_theta_true_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal
marginal_spread_theta_resp_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal
marginal_spread_resp_err_all          = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal
marginal_spread_confidence_report_all = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal

marginal_spread_vals = groupValues(sidx);

for i=1:n_uncertainty_levels
        
    grpIdx = sidx(i);
    fltIDx = ( G == grpIdx) ;
    
    % WRT actual stim ori
    grpRawErr_Flt   = filteredTable.rawOriError(fltIDx);
    grpStimOri      = filteredTable.stimOri(fltIDx);
    grpOriResp      = filteredTable.reportedOri(fltIDx);
    grpReportedConf = filteredTable.reportedConf(fltIDx);
    
    marginal_spread_theta_true_all(i, :)          = grpStimOri(:);
    marginal_spread_theta_resp_all(i, :)          = grpOriResp(:);
    marginal_spread_confidence_report_all(i, :)   = grpReportedConf(:);
    marginal_spread_resp_err_all(i, :)            = grpRawErr_Flt(:);
end

retData.marginal_spread_vals                   = marginal_spread_vals;
retData.marginal_spread_theta_true_all         = marginal_spread_theta_true_all;
retData.marginal_spread_theta_resp_all         = marginal_spread_theta_resp_all;
retData.marginal_spread_resp_err_all           = marginal_spread_resp_err_all;
retData.marginal_spread_confidence_report_all  = marginal_spread_confidence_report_all;


% Spread and duration
tableVals      = [data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur]; % Combine relevant columns from the table
rowsToKeep     = ismember(tableVals, combinations.sp_dur, 'rows');                % Find rows matching any combination in factorialCombinations1
filteredTable  = data.dat(rowsToKeep, :); % Filter table

% 1: contrast groups ==============================================
[G, groupValues]     = findgroups(filteredTable.stimDur);
grpIdxes             = unique(G);
n_uncertainty_levels = numel(grpIdxes);

sidx = getSortedGrpIdxes(n_uncertainty_levels, G, grpIdxes, filteredTable);

marginal_dur_theta_true_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both dur level no of trials must be equal
marginal_dur_theta_resp_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both dur level no of trials must be equal
marginal_dur_resp_err_all          = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both dur level no of trials must be equal
marginal_dur_confidence_report_all = zeros(n_uncertainty_levels, numel(filteredTable.stimContrast)/ 2); % Both both dur level no of trials must be equal

marginal_dur_vals = groupValues(sidx);

for i=1:n_uncertainty_levels
        
    grpIdx = sidx(i);
    fltIDx = ( G == grpIdx) ;
    
    % WRT actual stim ori
    grpRawErr_Flt   = filteredTable.rawOriError(fltIDx);
    grpStimOri      = filteredTable.stimOri(fltIDx);
    grpOriResp      = filteredTable.reportedOri(fltIDx);
    grpReportedConf = filteredTable.reportedConf(fltIDx);
    
    marginal_dur_theta_true_all(i, :)          = grpStimOri(:);
    marginal_dur_theta_resp_all(i, :)          = grpOriResp(:);
    marginal_dur_confidence_report_all(i, :)   = grpReportedConf(:);
    marginal_dur_resp_err_all(i, :)            = grpRawErr_Flt(:);
end

retData.marginal_dur_vals                   = marginal_dur_vals;
retData.marginal_dur_theta_true_all         = marginal_dur_theta_true_all;
retData.marginal_dur_theta_resp_all         = marginal_dur_theta_resp_all;
retData.marginal_dur_resp_err_all           = marginal_dur_resp_err_all;
retData.marginal_dur_confidence_report_all  = marginal_dur_confidence_report_all;


% 1: spread groups ==============================================
[G, groupValues]     = findgroups(filteredTable.stimSpread);
grpIdxes             = unique(G);
n_uncertainty_levels = numel(grpIdxes);

sidx = getSortedGrpIdxes(n_uncertainty_levels, G, grpIdxes, filteredTable);

marginal_spread_theta_true_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal
marginal_spread_theta_resp_all        = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal
marginal_spread_resp_err_all          = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal
marginal_spread_confidence_report_all = zeros(n_uncertainty_levels, numel(filteredTable.stimSpread)/ 2); % Both both spread level no of trials must be equal

marginal_spread_vals = groupValues(sidx);

for i=1:n_uncertainty_levels
        
    grpIdx = sidx(i);
    fltIDx = ( G == grpIdx) ;
    
    % WRT actual stim ori
    grpRawErr_Flt   = filteredTable.rawOriError(fltIDx);
    grpStimOri      = filteredTable.stimOri(fltIDx);
    grpOriResp      = filteredTable.reportedOri(fltIDx);
    grpReportedConf = filteredTable.reportedConf(fltIDx);
    
    marginal_spread_theta_true_all(i, :)          = grpStimOri(:);
    marginal_spread_theta_resp_all(i, :)          = grpOriResp(:);
    marginal_spread_confidence_report_all(i, :)   = grpReportedConf(:);
    marginal_spread_resp_err_all(i, :)            = grpRawErr_Flt(:);
end

retData.marginal_spread_vals2                   = marginal_spread_vals;
retData.marginal_spread_theta_true_all2         = marginal_spread_theta_true_all;
retData.marginal_spread_theta_resp_all2         = marginal_spread_theta_resp_all;
retData.marginal_spread_resp_err_all2           = marginal_spread_resp_err_all;
retData.marginal_spread_confidence_report_all2  = marginal_spread_confidence_report_all;



end



function sidx = getSortedGrpIdxes(n_uncertainty_levels, G, grpIdxes, filteredTable)

fitStds = zeros(1, n_uncertainty_levels);

for i=1:n_uncertainty_levels
    grpIdx          = grpIdxes(i);
    grpRawErr_Flt   = filteredTable.rawOriError(G == grpIdx);
    
    assert(all(~isnan(grpRawErr_Flt))) % there should be no nan
    pd = fitdist(grpRawErr_Flt, 'Normal'); % there should be non nan so need to filter nans

    sigma = pd.sigma;
    fitStds(i) = sigma;
end

[~, idx] = sort(fitStds);
sidx = grpIdxes(idx);

end
