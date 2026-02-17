clear all
close all

%% plot error histogram for each subject
figure
oriBins = -90:3:90;

for i = 1:20
    if i < 10
        data = load(sprintf('S0%d.mat', i));
    else
        data = load(sprintf('S%d.mat', i));
    end
    
    % respErr = data.RespErrDebias;
    stimOri = data.StimOri;
    respOri = data.RespOri;
    err     = respOri - stimOri;
    respErr = mod(err + 90, 180) - 90;
    
    subplot(4, 5, i)
    histogram(respErr, oriBins)
    title(sprintf('S%d.mat, STD: %.2f', i, std(respErr, 'omitnan') ))
    
end


%%
% I should probably not be using debaised error and any filtered error

data = load("S12.mat");
stimOri = data.StimOri;
respOri = data.RespOri;
% respErr = data.RespErrDebias;
conf    = data.Conf; % why does it go from -1 to 1

err     = respOri - stimOri;
respErr = mod(err + 90, 180) - 90;

mu = mean(respErr, 'omitnan');                    % Mean of the array
sigma = std(respErr, 'omitnan');                  % Standard deviation
% fltIdx = abs(respErr - mu) <= 3*sigma;

% stimOri = stimOri(fltIdx);
% respOri = respOri(fltIdx);
% conf    = conf(fltIdx);
% respErr = respErr(fltIdx); 

% respErr = (respOri - stimOri);
% respErr = mod(respErr + 90, 180) - 90; % Find minimum acute angle error

bins = 0:18:180;
binsCenters = ( bins(1:end-1) + bins(2:end) ) / 2;

% 1. Assign each stimOri to a bin index
[~, ~, binIdx] = histcounts(stimOri, bins);

% 2. Group respOri according to those bin indices
numBins = length(bins)-1;
groupedRespErr = cell(numBins,1); 
groupedStim    = cell(numBins,1);
groupedConf    = cell(numBins,1);
groupStds      = zeros(numBins,1);
groupStdsHC    = zeros(numBins,1);
groupStdsLC    = zeros(numBins,1);

for b = 1:numBins
    inBin = binIdx == b;
    groupedRespErr{b} = respErr(inBin);
    groupedStim{b}    = stimOri(inBin);
    groupedConf       = conf(inBin);
    
    groupStds(b) = std(groupedRespErr{b}, 'omitnan');
    data = groupedRespErr{b};
    groupStdsHC(b) = std(data(groupedConf > 0), 'omitnan');
    groupStdsLC(b) = std(data(groupedConf < 0), 'omitnan');
end


[~, idx] = sort(groupStds);

figure
subplot(2, 2, 1)
bar(1:numel(binsCenters), groupStds(idx)) % binsCenters(idx)

subplot(2, 2, 2)
hold on
plot(1:numel(binsCenters), groupStdsHC(idx), DisplayName="HC")
plot(1:numel(binsCenters), groupStdsLC(idx), DisplayName="LC")
hold off
legend


% figure
% for b = 1:numBins
%     data = groupedRespErr{b};
% 
%     subplot(2, 5, b)
%     histogram(data)
% end

% Fit parameters and do simulation in these value range (can you see the multiplicative effect)
