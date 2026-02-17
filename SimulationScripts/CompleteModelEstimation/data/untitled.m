m1 = load('ModelFitMetrics2.mat');
m2 = load('ModelFitMetrics3.mat');

m3 = [m1.dataToSave.fitInfo; m2.dataToSave.fitInfo];
dataToSave.groundTruth = m1.dataToSave.groundTruth;
dataToSave.fitInfo = m3;

save("ModelFitMetrics.mat", "dataToSave");