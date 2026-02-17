close all
clear all

data = load("ExpDataModelComparison.mat");
fValsCov = data.data.fValsCov;
fValsInd = data.data.fValsInd;

BinEdges = 0:10:1000;

figure

subplot(2, 2, 1)
hold on
histogram(fValsCov, BinEdges=BinEdges, DisplayName="Cov model") %BinEdges=BinEdges,
histogram(fValsInd, BinEdges=BinEdges, DisplayName="Independent model") %BinEdges=BinEdges,
xlabel("fvals")
ylabel("count")
title("Exp data (Raw data)")
legend
hold off

data = load("FltExpDataModelComparison.mat");
fValsCov = data.data.fValsCov;
fValsInd = data.data.fValsInd;

BinEdges = 0:1:100;

subplot(2, 2, 2)
hold on
histogram(fValsCov, BinEdges=BinEdges, DisplayName="Cov model") %BinEdges=BinEdges,
histogram(fValsInd, BinEdges=BinEdges, DisplayName="Independent model") %BinEdges=BinEdges,
% histogram(fValsCov, DisplayName="Cov model") %BinEdges=BinEdges,
% histogram(fValsInd, DisplayName="Independent model") %BinEdges=BinEdges,
xlabel("fvals")
ylabel("count")
title("Exp data (flt data)")
legend
hold off
