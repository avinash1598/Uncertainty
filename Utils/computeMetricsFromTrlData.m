function metrics = computeMetricsFromTrlData(trlErrors, trlConfReports, trlUncertaintyLevels)

uniqLevels = sort(unique(trlUncertaintyLevels));

mads    = zeros(1, numel(uniqLevels));
mads_HC = zeros(1, numel(uniqLevels));
mads_LC = zeros(1, numel(uniqLevels));

for i = 1:numel(uniqLevels)
    fltIdx  = trlUncertaintyLevels == uniqLevels(i);
    mads(i) = mad(trlErrors(fltIdx), 1); % wrt median

    % HC
    fltIdx     = ( trlUncertaintyLevels == uniqLevels(i) ) & ( trlConfReports == 1 );
    mads_HC(i) = mad(trlErrors(fltIdx), 1); % wrt median
    
    % LC
    fltIdx     = ( trlUncertaintyLevels == uniqLevels(i) ) & ( trlConfReports == 0 );
    mads_LC(i) = mad(trlErrors(fltIdx), 1); % wrt median

end

metrics.mads    = mads;
metrics.mads_HC = mads_HC;
metrics.mads_LC = mads_LC;

end
