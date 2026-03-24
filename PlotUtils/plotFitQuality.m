function plotFitQuality(data, fitData)

set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',1, ...
    'defaultAxesBox','off', ...
    'defaultAxesTickDir','out', ...
    'defaultLineLineWidth',1.5, ...
    'defaultAxesLabelFontSizeMultiplier',1.1, ...
    'defaultAxesTitleFontWeight','normal', ...
    'defaultLegendBox','off');

titleTxt = ["Cov Model (Reduced)", ...
            "Cov Model (Full)", ...
            "Cov Model (Full Jumbo)", ...
            "Ind Model (Reduced)", ...
            "Ind Model (Full)", ...
            "Ind Model (Full Jumbo)"];

keys = ["resultReducedCov", ...
        "resultFullCov", ...
        "resultFullJumboCov", ...
        "resultReducedInd", ...
        "resultFullInd", ...
        "resultFullJumboInd"];

keysSim = ["resultReducedCovSim", ...
           "resultFullCovSim", ...
           "resultFullJumboCovSim", ...
           "resultReducedIndSim", ...
           "resultFullIndSim", ...
           "resultFullJumboIndSim"];

% BinEdges = 11.1:0.01:11.5;
BinEdges = 7.1e4:20:7.3e4;

% IS the fit quality bad (6 plots)?
% Plot cov (3 plots), Ind (3 plots)
figure

for i=1:numel(keys)
    %fvals    = log( fitData.(keys(i)).f );
    %fvalsSim = log( fitData.(keysSim(i)).f );
    fvals    = ( fitData.(keys(i)).f );
    fvalsSim = ( fitData.(keysSim(i)).f );
    
    top100 = sort(fvals, 'ascend'); top100 = top100(1:200); fvals = top100;
    top100 = sort(fvalsSim, 'ascend'); top100 = top100(1:200); fvalsSim = top100;
    
    % fvals = fvals / numel( data.theta_true_all );
    %BinEdges_ = 7e4:50:7.5e4;

    subplot(2, 3, i)
    hold on
    histogram(fvals,  BinEdges, DisplayName="Data"); %BinEdges,
    histogram(fvalsSim,  BinEdges, DisplayName="Simulation");
    %xlabel("fvals (log)")
    xlabel("fvals")
    ylabel("count")
    title(titleTxt(i))
    legend
    hold off

    % p = ranksum(fvals, fvalsSim)
end

set(gca, 'XScale', 'log');

% Does complexity of the model helps with the fitting quality? (cov, ind - 4 plots - 2 real data, 2 simulation)
figure

subplot(2, 2, 1)
hold on
histogram((fitData.resultReducedCov.f), BinEdges, DisplayName="Reduced");
histogram((fitData.resultFullCov.f), BinEdges, DisplayName="Full");
histogram((fitData.resultFullJumboCov.f), BinEdges, DisplayName="Full Jumbo");
xlabel("fvals log")
ylabel("count")
title("Cov model: Fit quality (Data)")
legend
hold off

subplot(2, 2, 2)
hold on
histogram((fitData.resultReducedInd.f), BinEdges, DisplayName="Reduced");
histogram((fitData.resultFullInd.f), BinEdges, DisplayName="Full");
histogram((fitData.resultFullJumboInd.f), BinEdges, DisplayName="Full Jumbo");
xlabel("fvals log")
ylabel("count")
title("Ind model: Fit quality (Data)")
legend
hold off

subplot(2, 2, 3)
hold on
histogram((fitData.resultReducedCovSim.f), BinEdges, DisplayName="Reduced");
histogram((fitData.resultFullCovSim.f), BinEdges, DisplayName="Full");
histogram((fitData.resultFullJumboCovSim.f), BinEdges, DisplayName="Full Jumbo");
xlabel("fvals log")
ylabel("count")
title("Cov model: Fit quality (Simulation)")
legend
hold off

subplot(2, 2, 4)
hold on
histogram((fitData.resultReducedIndSim.f), BinEdges, DisplayName="Reduced");
histogram((fitData.resultFullIndSim.f), BinEdges, DisplayName="Full");
histogram((fitData.resultFullJumboIndSim.f), BinEdges, DisplayName="Full Jumbo");
xlabel("fvals log")
ylabel("count")
title("Ind model: Fit quality (Simulation)")
legend
hold off

set(gca, 'XScale', 'log');

% Is one model better than the other? 3 plots
figure

subplot(2, 3, 1)
hold on
histogram((fitData.resultReducedCov.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultReducedInd.f), BinEdges, DisplayName="Ind");
xlabel("fvals log")
ylabel("count")
title("Reduced Model (fit to data)")
legend
hold off

subplot(2, 3, 2)
hold on
histogram((fitData.resultFullCov.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultFullInd.f), BinEdges, DisplayName="Ind");
xlabel("fvals log")
ylabel("count")
title("Full Model (fit to data)")
legend
hold off

subplot(2, 3, 3)
hold on
histogram((fitData.resultFullJumboCov.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultFullJumboInd.f), BinEdges, DisplayName="Ind");
xlabel("fvals log")
ylabel("count")
title("Full Jumbo Model (fit to data)")
legend
hold off

subplot(2, 3, 4)
hold on
histogram((fitData.resultReducedCovSim.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultReducedIndSim.f), BinEdges, DisplayName="Ind");
xlabel("fvals log")
ylabel("count")
title("Reduced Model (simulated data)")
legend
hold off

subplot(2, 3, 5)
hold on
histogram((fitData.resultFullCovSim.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultFullIndSim.f), BinEdges, DisplayName="Ind");
xlabel("fvals log")
ylabel("count")
title("Full Model (simulated data)")
legend
hold off

subplot(2, 3, 6)
hold on
histogram((fitData.resultFullJumboCovSim.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultFullJumboIndSim.f), BinEdges, DisplayName="Ind");
xlabel("fvals log")
ylabel("count")
title("Full Jumbo Model (simulated data)")
legend
hold off

set(gca, 'XScale', 'log');

end
