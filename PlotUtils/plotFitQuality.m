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
% BinEdges = 7.1e4:20:7.17e4;
%BinEdges = 7.10e4:20:7.20e4;
BinEdges = 4.00e4:20:7.00e4;

% IS the fit quality bad (6 plots)?
% Plot cov (3 plots), Ind (3 plots)
figure

for i=1:numel(keys)
    %fvals    = log( fitData.(keys(i)).f );
    %fvalsSim = log( fitData.(keysSim(i)).f );
    fvals    = ( fitData.(keys(i)).f );
    fvalsSim = ( fitData.(keysSim(i)).f );
    
    %top100 = sort(fvals, 'ascend'); top100 = top100(1:200); fvals = top100;
    %top100 = sort(fvalsSim, 'ascend'); top100 = top100(1:200); fvalsSim = top100;
    
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
    %xticks(linspace(min([fvals fvalsSim]), max([fvals fvalsSim]), 4)); 

    % p = ranksum(fvals, fvalsSim)
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.XTickMode = 'auto';
    ax.XTickLabelMode = 'auto';

end

% set(gca, 'XScale', 'log');

% Does complexity of the model helps with the fitting quality? (cov, ind - 4 plots - 2 real data, 2 simulation)
figure

subplot(2, 2, 1)
hold on
histogram((fitData.resultReducedCov.f), BinEdges, DisplayName="Reduced");
histogram((fitData.resultFullCov.f), BinEdges, DisplayName="Full");
histogram((fitData.resultFullJumboCov.f), BinEdges, DisplayName="Full Jumbo");
xlabel("fvals")
ylabel("count")
title("Cov model: Fit quality (Data)")
legend
hold off

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

subplot(2, 2, 2)
hold on
histogram((fitData.resultReducedInd.f), BinEdges, DisplayName="Reduced");
histogram((fitData.resultFullInd.f), BinEdges, DisplayName="Full");
histogram((fitData.resultFullJumboInd.f), BinEdges, DisplayName="Full Jumbo");
xlabel("fvals")
ylabel("count")
title("Ind model: Fit quality (Data)")
legend
hold off

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

BinEdges2 = -600:10:600;

n = 1728;        % number of trials
k_red = 10;      % parameters (Cov model)
k_full = 12;      % parameters (Ind model)
nll_red = fitData.resultReducedCov.f;
nll_full = fitData.resultFullCov.f;
[delAIC, delBIC] = computeAIC_BIC(nll_full, nll_red, k_full, k_red, n);

subplot(2,2,3); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (Full - Reduced)")
title('Cov Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

n = 1728;        % number of trials
k_red = 11;      % parameters (Cov model)
k_full = 13;      % parameters (Ind model)
nll_red = fitData.resultReducedInd.f;
nll_full = fitData.resultFullInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_full, nll_red, k_full, k_red, n);

subplot(2,2,4); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (Full - Reduced)")
title('Ind Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

% subplot(2, 2, 3)
% hold on
% histogram((fitData.resultReducedCovSim.f), BinEdges, DisplayName="Reduced");
% histogram((fitData.resultFullCovSim.f), BinEdges, DisplayName="Full");
% histogram((fitData.resultFullJumboCovSim.f), BinEdges, DisplayName="Full Jumbo");
% xlabel("fvals")
% ylabel("count")
% title("Cov model: Fit quality (Simulation)")
% legend
% hold off
% 
% subplot(2, 2, 4)
% hold on
% histogram((fitData.resultReducedIndSim.f), BinEdges, DisplayName="Reduced");
% histogram((fitData.resultFullIndSim.f), BinEdges, DisplayName="Full");
% histogram((fitData.resultFullJumboIndSim.f), BinEdges, DisplayName="Full Jumbo");
% xlabel("fvals")
% ylabel("count")
% title("Ind model: Fit quality (Simulation)")
% legend
% hold off

% set(gca, 'XScale', 'log');

% Is one model better than the other? 3 plots
figure

subplot(3, 3, 1)
hold on
histogram((fitData.resultReducedCov.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultReducedInd.f), BinEdges, DisplayName="Ind");
xlabel("fvals")
ylabel("count")
title("Reduced Model (fit to data)")
legend
hold off

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';


subplot(3, 3, 2)
hold on
histogram((fitData.resultFullCov.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultFullInd.f), BinEdges, DisplayName="Ind");
xlabel("fvals")
ylabel("count")
title("Full Model (fit to data)")
legend
hold off

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

subplot(3, 3, 3)
hold on
histogram((fitData.resultFullJumboCov.f), BinEdges, DisplayName="Cov");
histogram((fitData.resultFullJumboInd.f), BinEdges, DisplayName="Ind");
xlabel("fvals")
ylabel("count")
title("Full Jumbo Model (fit to data)")
legend
hold off

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

BinEdges2 = -300:10:300;

n = 1728;        % number of trials
k_cov = 10;      % parameters (Cov model)
k_ind = 11;      % parameters (Ind model)
nll_cov = fitData.resultReducedCov.f;
nll_ind = fitData.resultReducedInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,4); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
% histogram(delBIC, 'DisplayName','delta BIC');
xlabel("delta (COV - IND)")
title('Reduced Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

n = 1728;        % number of trials
k_cov = 12;      % parameters (Cov model)
k_ind = 13;      % parameters (Ind model)
nll_cov = fitData.resultFullCov.f;
nll_ind = fitData.resultFullInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,5); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
% histogram(delBIC, 'DisplayName','delta BIC');
xlabel("delta (COV - IND)")
title('Full Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

n = 1728;        % number of trials
k_cov = 17;      % parameters (Cov model)
k_ind = 18;      % parameters (Ind model)
nll_cov = fitData.resultFullJumboCov.f;
nll_ind = fitData.resultFullJumboInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,6); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
% histogram(delBIC, 'DisplayName','delta BIC');
xlabel("delta (COV - IND)")
title('Full Jumbo Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

n = 1728;        % number of trials
k_cov = 10;      % parameters (Cov model)
k_ind = 11;      % parameters (Ind model)
nll_cov = fitData.resultReducedCov.f;
nll_ind = fitData.resultReducedInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,7); hold on;
% histogram(delAIC, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (COV - IND)")
title('Reduced Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

n = 1728;        % number of trials
k_cov = 12;      % parameters (Cov model)
k_ind = 13;      % parameters (Ind model)
nll_cov = fitData.resultFullCov.f;
nll_ind = fitData.resultFullInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,8); hold on;
% histogram(delAIC, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (COV - IND)")
title('Full Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

n = 1728;        % number of trials
k_cov = 17;      % parameters (Cov model)
k_ind = 18;      % parameters (Ind model)
nll_cov = fitData.resultFullJumboCov.f;
nll_ind = fitData.resultFullJumboInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,9); hold on;
% histogram(delAIC, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (COV - IND)")
title('Full Jumbo Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';


% subplot(2, 3, 4)
% hold on
% histogram((fitData.resultReducedCovSim.f), BinEdges, DisplayName="Cov");
% histogram((fitData.resultReducedIndSim.f), BinEdges, DisplayName="Ind");
% xlabel("fvals")
% ylabel("count")
% title("Reduced Model (simulated data)")
% legend
% hold off
% 
% subplot(2, 3, 5)
% hold on
% histogram((fitData.resultFullCovSim.f), BinEdges, DisplayName="Cov");
% histogram((fitData.resultFullIndSim.f), BinEdges, DisplayName="Ind");
% xlabel("fvals")
% ylabel("count")
% title("Full Model (simulated data)")
% legend
% hold off
% 
% subplot(2, 3, 6)
% hold on
% histogram((fitData.resultFullJumboCovSim.f), BinEdges, DisplayName="Cov");
% histogram((fitData.resultFullJumboIndSim.f), BinEdges, DisplayName="Ind");
% xlabel("fvals")
% ylabel("count")
% title("Full Jumbo Model (simulated data)")
% legend
% hold off

% set(gca, 'XScale', 'log');

end

function [delAIC, delBIC] = computeAIC_BIC(nll1, nll2, k1, k2, n)

% AIC
AIC_1 = 2*k1 + 2*nll1;
AIC_2 = 2*k2 + 2*nll2;

% BIC
BIC_1 = k1*log(n) + 2*nll1;
BIC_2 = k2*log(n) + 2*nll2;

delAIC = AIC_1 - AIC_2;
delBIC = BIC_1 - BIC_2;

end
