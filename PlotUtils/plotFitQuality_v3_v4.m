function plotFitQuality_v3_v4(data, fitData1, fitData2, fitData3)

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
BinEdges = 7.105e4:2:7.12e4; % Akash
% BinEdges = 7.09e4:2:7.14e4; % Jonathan
% BinEdges = 7.1e4:2:7.11e4; % Tien
% BinEdges = 7.08e4:2:7.113e4; % Yichao

% IS the fit quality bad (6 plots)?
% Plot cov (3 plots), Ind (3 plots)

% set(gca, 'XScale', 'log');

% Is one model better than the other? 3 plots
figure

subplot(3, 3, 1)
hold on
histogram((fitData1.resultReducedCov.f), BinEdges, DisplayName="Cov");
histogram((fitData1.resultReducedInd.f), BinEdges, DisplayName="Ind");
histogram((fitData2.result.f), BinEdges, DisplayName="ADD+MULT");
histogram((fitData3.result.f), BinEdges, DisplayName="ADD+MULT (FULL)");
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
histogram((fitData1.resultFullCov.f), BinEdges, DisplayName="Cov");
histogram((fitData1.resultFullInd.f), BinEdges, DisplayName="Ind");
histogram((fitData2.result.f), BinEdges, DisplayName="ADD+MULT");
histogram((fitData3.result.f), BinEdges, DisplayName="ADD+MULT (FULL)");
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
histogram((fitData1.resultFullJumboCov.f), BinEdges, DisplayName="Cov");
histogram((fitData1.resultFullJumboInd.f), BinEdges, DisplayName="Ind");
histogram((fitData2.result.f), BinEdges, DisplayName="ADD+MULT");
histogram((fitData3.result.f), BinEdges, DisplayName="ADD+MULT (FULL)");
xlabel("fvals")
ylabel("count")
title("Full Jumbo Model (fit to data)")
legend
hold off

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

% BinEdges2 = -300:10:300; % Akash
BinEdges2 = -600:10:600;

n = 1728;        % number of trials
k_cov = 10;      % parameters (Cov model)
k_ind = 11;      % parameters (Ind model)
nll_cov = fitData1.resultReducedCov.f;
nll_ind = fitData1.resultReducedInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,4); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
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
nll_cov = fitData1.resultFullCov.f;
nll_ind = fitData1.resultFullInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,5); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
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
nll_cov = fitData1.resultFullJumboCov.f;
nll_ind = fitData1.resultFullJumboInd.f;
[delAIC, delBIC] = computeAIC_BIC(nll_cov, nll_ind, k_cov, k_ind, n);

subplot(3,3,6); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
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

%% AIC BIC for ADD+MULT Model

f1 = fitData1.resultReducedCov.f;
f2 = fitData3.result.f;
top100 = sort(f1, 'ascend'); top100 = top100(1:100); f1 = top100;
top100 = sort(f2, 'ascend'); top100 = top100(1:100); f2 = top100;
    
n = 1728;        % number of trials
k_cov = 10;      % parameters (Cov model)
k_add_mult = 18; %16;      % parameters (Ind model)
nll_cov = f1;
nll_add_mult = f2; % fitData2.result.f;
[delAIC, delBIC] = computeAIC_BIC(nll_add_mult, nll_cov, k_add_mult, k_cov, n);

subplot(3,3,7); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (ADD+MULT - COV)")
title('Reduced Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';

f1 = fitData1.resultFullCov.f;
f2 = fitData3.result.f;
top100 = sort(f1, 'ascend'); top100 = top100(1:100); f1 = top100;
top100 = sort(f2, 'ascend'); top100 = top100(1:100); f2 = top100;
 
n = 1728;        % number of trials
k_cov = 12;      % parameters (Cov model)
k_add_mult = 18;  %16    % parameters (Ind model)
nll_cov = f1;
nll_add_mult = f2; % fitData2.result.f;
[delAIC, delBIC] = computeAIC_BIC(nll_add_mult, nll_cov, k_add_mult, k_cov, n);

subplot(3,3,8); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (ADD+MULT - COV)")
title('Full Model');
legend;
xline(0, LineStyle="--")
% xlim([-300, 300])

ax = gca;
ax.XAxis.Exponent = 0;
ax.XTickMode = 'auto';
ax.XTickLabelMode = 'auto';


f1 = fitData1.resultFullJumboCov.f;
f2 = fitData3.result.f;
top100 = sort(f1, 'ascend'); top100 = top100(1:100); f1 = top100;
top100 = sort(f2, 'ascend'); top100 = top100(1:100); f2 = top100;
assert(size(fitData1.resultFullJumboCov.x, 2) == 17)
assert(size(fitData3.result.x, 2) == 18)
n = 1728;        % number of trials
k_cov = 17;      % parameters (Cov model)
k_add_mult = 18; %16      % parameters (Ind model)
nll_cov = f1;
nll_add_mult = f2; %fitData2.result.f;
[delAIC, delBIC] = computeAIC_BIC(nll_add_mult, nll_cov, k_add_mult, k_cov, n);

subplot(3,3,9); hold on;
histogram(delAIC, BinEdges2, 'DisplayName',' delta AIC');
histogram(delBIC, BinEdges2, 'DisplayName','delta BIC');
xlabel("delta (ADD+MULT - COV)")
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
% histogram((fitData1.resultReducedCovSim.f), BinEdges, DisplayName="Cov");
% histogram((fitData1.resultReducedIndSim.f), BinEdges, DisplayName="Ind");
% xlabel("fvals")
% ylabel("count")
% title("Reduced Model (simulated data)")
% legend
% hold off
% 
% subplot(2, 3, 5)
% hold on
% histogram((fitData1.resultFullCovSim.f), BinEdges, DisplayName="Cov");
% histogram((fitData1.resultFullIndSim.f), BinEdges, DisplayName="Ind");
% xlabel("fvals")
% ylabel("count")
% title("Full Model (simulated data)")
% legend
% hold off
% 
% subplot(2, 3, 6)
% hold on
% histogram((fitData1.resultFullJumboCovSim.f), BinEdges, DisplayName="Cov");
% histogram((fitData1.resultFullJumboIndSim.f), BinEdges, DisplayName="Ind");
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
