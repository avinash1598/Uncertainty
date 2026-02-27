function binnedData = buildBinnedData( ...
    n_levels, ...
    errBins, ...
    trlErrs, ...
    trlConfReports, ...
    trlUncertaintyLevels, ...
    trlStimOris, ...
    orientations, ...
    binnedByOri)

if nargin < 8 && isempty(binnedByOri)
    binnedByOri = false;
end

if binnedByOri
    binnedData = GetBinnedDataByOri( ...
        n_levels, ...
        errBins, ...
        trlErrs, ...
        trlConfReports, ...
        trlUncertaintyLevels, ...
        trlStimOris, ...
        orientations);
else
    binnedData = GetBinnedData( ...
        n_levels, ...
        errBins, ...
        trlErrs, ...
        trlConfReports, ...
        trlUncertaintyLevels);
end

end

function binnedData = GetBinnedData(n_levels, errBins, trlErrs, ...
    trlConfReports, trlUncertaintyLevels)

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

function binnedData = GetBinnedDataByOri(n_levels, errBins, trlErrs, ...
    trlConfReports, trlUncertaintyLevels, trlStimOris, orientations)

% Get PDFs from data for HC and LC
binned_err_LC = zeros( n_levels, numel(orientations), numel(errBins) );
binned_err_HC = zeros( n_levels, numel(orientations), numel(errBins) );

for i=1:n_levels
    for j = 1:numel(orientations)
        ori = orientations(j);
        fltIdx = (trlUncertaintyLevels == i) & (trlStimOris == ori);
        
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
        
        binned_err_LC(i, j, :) = countLC;
        binned_err_HC(i, j, :) = countHC;
    end
end

binnedData.binned_err_LC = binned_err_LC;
binnedData.binned_err_HC = binned_err_HC;

end



% function binnedData = buildBinnedData(n_levels, errBins, trlErrs, trlConfReports, trlUncertaintyLevels)
% 
% % Get PDFs from data for HC and LC
% binned_err_LC = zeros( n_levels, numel(errBins) );
% binned_err_HC = zeros( n_levels, numel(errBins) );
% 
% for i=1:n_levels
%     
%     fltIdx = trlUncertaintyLevels == i;
% 
%     cR = trlConfReports(fltIdx);
%     fltErr = trlErrs(fltIdx);
% 
%     dataHC = fltErr(cR == 1);
%     dataLC = fltErr(cR == 0);
%     
%     centers = errBins;
%     binWidth = mean(diff(centers));
%     edges = [centers - binWidth/2, centers(end) + binWidth/2];
%     
%     [countHC, ~] = histcounts(dataHC, ...
%         'Normalization', 'count', ...
%         'BinEdges', edges);
% 
%     [countLC, ~] = histcounts(dataLC, ...
%         'Normalization', 'count', ...
%         'BinEdges', edges);
%     
%     binned_err_LC(i, :) = countLC;
%     binned_err_HC(i, :) = countHC;
%     
% end
% 
% binnedData.binned_err_LC = binned_err_LC;
% binnedData.binned_err_HC = binned_err_HC;
% 
% end