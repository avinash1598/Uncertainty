function trlData = convertToTrialData(data)

grpOriErr            = data.resp_err_all; 
confReport           = data.confidence_report_all;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr    = reshape(grpOriErr, n_uncertainty_levels, []);
confReport   = reshape(confReport, n_uncertainty_levels, []);

levels               = 1:n_uncertainty_levels;
uncertainty_levels   = repmat(levels', [1 size(grpOriErr, 2)]);

% Compute MADs and Stds
y_std   = std(grpOriErr, 0, 2);
y_mad = mad(grpOriErr, 1, 2);

HC_idx = confReport == 1;
LC_idx = confReport == 0;

resp_HC = grpOriErr;
resp_HC(~HC_idx) = NaN;

resp_LC = grpOriErr;
resp_LC(~LC_idx) = NaN;

y_HC_std = std(resp_HC, 0, 2, 'omitnan');
y_HC_mad = mad(resp_HC, 1, 2);

y_LC_std = std(resp_LC, 0, 2, 'omitnan');
y_LC_mad = mad(resp_LC, 1, 2);


% Faltten data structures
trlErrors            = grpOriErr(:);
trlConfReports       = confReport(:);
trlUncertaintyLevels = uncertainty_levels(:);

trlData.grpOriErr            = grpOriErr;
trlData.n_uncertainty_levels = n_uncertainty_levels;
trlData.trlErrors            = trlErrors;
trlData.trlConfReports       = trlConfReports;
trlData.trlUncertaintyLevels = trlUncertaintyLevels;

trlData.y_std     = y_std;
trlData.y_HC_std  = y_HC_std;
trlData.y_LC_std  = y_LC_std;

trlData.y_mad     = y_mad;
trlData.y_HC_mad  = y_HC_mad;
trlData.y_LC_mad  = y_LC_mad;

end
