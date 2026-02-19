function binnedData = buildBinnedData(n_levels, errBins, trlErrs, trlConfReports, trlUncertaintyLevels)

% Get PDFs from data for HC and LC
binned_err_LC = zeros( n_levels, numel(errBins) );
binned_err_HC = zeros( n_levels, numel(errBins) );

for i=1:n_levels
    
    fltIdx = trlUncertaintyLevels == i;

    cR = trlConfReports(fltIdx);
    fltErr = trlErrs(fltIdx);

    dataHC = fltErr(cR == 1);
    dataLC = fltErr(cR == 0);
    
    centers = errBins;
    binWidth = mean(diff(centers));
    edges = [centers - binWidth/2, centers(end) + binWidth/2];
    
    [countHC, ~] = histcounts(dataHC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);

    [countLC, ~] = histcounts(dataLC, ...
        'Normalization', 'count', ...
        'BinEdges', edges);
    
    binned_err_LC(i, :) = countLC;
    binned_err_HC(i, :) = countHC;
    
end

binnedData.binned_err_LC = binned_err_LC;
binnedData.binned_err_HC = binned_err_HC;

end