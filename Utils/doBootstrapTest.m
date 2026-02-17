function [retData] = doBootstrapTest(groupA_H1, groupB_H1, groupA_H5, groupB_H5, noisecorr)
    % Bootstrap test to test if difference between group A and B at
    % hyperparameter 1 is greater than the difference at hypperparameter 2.
    
    % Set default value if noisecorr is not provided
    if nargin < 5 || isempty(noisecorr)
        noisecorr = false;  % or true, or any default value you want
    end

    % Bootstrap settings
    nBoot = 10000;
    diffDiff = zeros(nBoot, 1);
    
    % Bootstrap resampling loop
    for i = 1:nBoot
        % Resample with replacement
        A1 = groupA_H1(randi(numel(groupA_H1), numel(groupA_H1), 1));
        B1 = groupB_H1(randi(numel(groupB_H1), numel(groupB_H1), 1));
        A5 = groupA_H5(randi(numel(groupA_H5), numel(groupA_H5), 1));
        B5 = groupB_H5(randi(numel(groupB_H5), numel(groupB_H5), 1));
        
        % We only need to worry if absolute difference decreases over
        % hyperparameter values
        delta1 = abs(mean(A1) - mean(B1));
        delta5 = abs(mean(A5) - mean(B5));
        diffDiff(i) = delta1 - delta5;
    end
    
    % Observed difference of differences
    obsDelta1 = abs(mean(groupA_H1) - mean(groupB_H1));
    obsDelta5 = abs(mean(groupA_H5) - mean(groupB_H5));
    obsDiff = obsDelta1 - obsDelta5;
    
    % One-sided p-value: is Δ1 − Δ5 > 0?
    % pval is zero when the difference is greater than zero for all the
    % bootstrap iterations
    pval = mean(diffDiff <= 0);
    
    % If the simulation is using noise correlation then the test is to
    % check if difference is less than zero for all the iterations
    if noisecorr
        pval = mean(diffDiff >= 0);
    end
    
    retData.obsDiff = obsDiff;
    retData.pval = pval;
    
%     if pval > 0.05
%         fprintf("Bootstrap test : %f \n", pval)
%         fprintf("Observed diff : %f \n", obsDiff)
%     end
end