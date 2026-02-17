function trlData = convertToTrialData(data)

grpOriErr            = data.resp_err_all; 
confReport           = data.confidence_report_all;
n_uncertainty_levels = size(grpOriErr, 1);

grpOriErr    = reshape(grpOriErr, n_uncertainty_levels, []);
confReport   = reshape(confReport, n_uncertainty_levels, []);

levels               = 1:n_uncertainty_levels;
uncertainty_levels   = repmat(levels', [1 size(grpOriErr, 2)]);

% Faltten data structures
trlErrors            = grpOriErr(:);
trlConfReports       = confReport(:);
trlUncertaintyLevels = uncertainty_levels(:);

trlData.grpOriErr            = grpOriErr;
trlData.n_uncertainty_levels = n_uncertainty_levels;
trlData.trlErrors            = trlErrors;
trlData.trlConfReports       = trlConfReports;
trlData.trlUncertaintyLevels = trlUncertaintyLevels;

end
