close all
clear all

data = load("CovDataModelComparison_t1000.mat");
fValsCov = data.data.fValsCov;
fValsInd = data.data.fValsInd;

BinEdges = 0:0.1:2;

figure

subplot(2, 2, 1)
hold on
histogram(fValsCov, BinEdges=BinEdges, DisplayName="Cov model")
histogram(fValsInd,  BinEdges=BinEdges, DisplayName="Independent model")
xlabel("fvals")
ylabel("count")
title("Simulated data: cov (80000 trials)")
legend
hold off


data = load("IndDataModelComparison_t1000.mat");
fValsCov = data.data.fValsCov;
fValsInd = data.data.fValsInd;

BinEdges = 0:0.1:2;

subplot(2, 2, 2)
hold on
histogram(fValsCov, BinEdges=BinEdges, DisplayName="Cov model")
histogram(fValsInd,  BinEdges=BinEdges, DisplayName="Independent model")
xlabel("fvals")
ylabel("count")
title("Simulated data: independent (80000 trials)")
legend
hold off


data = load("CovDataModelComparison_t25.mat");
fValsCov = data.data.fValsCov;
fValsInd = data.data.fValsInd;

BinEdges = 9:0.1:11;

subplot(2, 2, 3)
hold on
histogram(fValsCov, BinEdges=BinEdges, DisplayName="Cov model")
histogram(fValsInd,  BinEdges=BinEdges, DisplayName="Independent model")
xlabel("fvals")
ylabel("count")
title("Simulated data: cov (2000 trials)")
legend
hold off


data = load("IndDataModelComparison_t25.mat");
fValsCov = data.data.fValsCov;
fValsInd = data.data.fValsInd;

BinEdges = 26.5:0.1:28;

subplot(2, 2, 4)
hold on
histogram(fValsCov, BinEdges=BinEdges, DisplayName="Cov model")
histogram(fValsInd,  BinEdges=BinEdges, DisplayName="Independent model")
xlabel("fvals")
ylabel("count")
title("Simulated data: Ind (2000 trials)")
legend
hold off

