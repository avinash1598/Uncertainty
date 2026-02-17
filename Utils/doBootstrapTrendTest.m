function [retData] = doBootstrapTrendTest(contrastLevels, groupHC, groupLC, noisecorr)
    % Boottrap test to check if the difference between groupA and groupB
    % decreases with increasing hyperparameter value
    % Pass groupA high confidence, groupB low confidence
    % Set default value if noisecorr is not provided
    if nargin < 4 || isempty(noisecorr)
        noisecorr = false;  % or true, or any default value you want
    end

    k = numel(contrastLevels); % No of contrast levels
    nBoot = 10000;
    bootstrapCorrelations = zeros(nBoot, 1);
    meanDiffs = zeros(1, k);
    
    assert((numel(contrastLevels) == size(groupHC, 1)) && (numel(contrastLevels) == size(groupLC, 1)) )

    % Step 1: Compute observed group differences
    meanDiffs(:) = abs( mean(groupHC, 2, 'omitnan') - mean(groupLC, 2, 'omitnan') );
    
    % Step 2: Difference is probably not linear so don't use regression
    % test. Use correlation instead.
    validIdx = ~isnan(contrastLevels) & ~isnan(meanDiffs);
    x1 = contrastLevels(validIdx);
    y1 = meanDiffs(validIdx);
    corr_obs = corr(x1(:), y1(:), 'Type', 'Spearman');
    
    % Step 3: Bootstrap
    for b = 1:nBoot
        bootDiffs = zeros(1, k);
        for i = 1:k
            x_ = groupHC(i, :);
            y_ = groupLC(i, :);
            x_ = x_(~isnan(x_));
            y_ = y_(~isnan(y_));
            
            % If any of the vector are empty then ignore this contrast
            % level for testing
            if ( numel(x_) == 0 ) || ( numel(y_) == 0 )
                bootDiffs(i) = NaN;
                continue
            end

            x_star = x_(randi(numel(x_), numel(x_), 1));
            y_star = y_(randi(numel(y_), numel(y_), 1));
            
            diff = mean(x_star) - mean(y_star);
            bootDiffs(i) = abs(diff);
        end

        x = contrastLevels(:);
        y = bootDiffs(:);
        valid = ~isnan(x) & ~isnan(y);

        bootstrapCorrelations(b) = corr(x(valid), y(valid), 'Type', 'Spearman');
    end
    
    % Step 4: Compute one-sided p-value (test for decreasing trend)
    % Perfect scenraio - pval=0: when all the correlations are negative
    pval = mean(bootstrapCorrelations >= 0);
    
    % If the simulation is using noise correlation then the test is to
    % check if all the correlations are positive (for Perfect scenraio - pval=0)
    if noisecorr
        pval = mean(bootstrapCorrelations <= 0);
    end
    
    retData.obs_corr = corr_obs;
    retData.pval = pval;

    if pval > 0.05
        fprintf("Bootstrap trend test : %f \n", pval)
        fprintf("Observed correlation : %f \n", corr_obs)
    end
end