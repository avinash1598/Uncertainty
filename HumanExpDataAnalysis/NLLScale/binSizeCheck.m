clc
clear all
close all

%Based on analysis it might be good idea to use binsize of 1 and multiply
%by binwidth to get to the true NLL

% 1) simulate data
N = 1e4; sigma = 50;
data = sigma * randn(N,1);
data = mod(data + 90, 180) - 90;

% 2) true NLL (analytic Gaussian)
L = 180; K = 5;
ks = -K:K;
pdf_wrap = @(x) sum(normpdf(x + ks*L, 0, sigma), 2);
trueNLL = -sum(log(pdf_wrap(data) + eps));
% trueNLL = -sum(log(normpdf(data,0,sigma)));

% KDE
[f, xi] = ksdensity(data, 'Function', 'pdf', 'Support', [-90 90]);
likelihood = interp1(xi, f, data, 'linear'); 
nllKDE = -log(likelihood);
nllKDE = sum(nllKDE(:));
fprintf("NLL KDE: %.4f True NLL: %.4f diff: %.4f\n", ...
    nllKDE, trueNLL, nllKDE - trueNLL)

figure
histogram(data, 100, Normalization="pdf")
hold on
plot(xi, f, DisplayName="fit", LineWidth=1.5)
hold off

% 3) binned NLL for different bin widths
% NLL artifically decreases because of binSize term in pdf
bws = [0.001 1 2 5 10 50];
for bw = bws
    edges = -90:bw:90; edges(end+1)=max(edges(end),90);
    p = histcounts(data, edges, 'Normalization','pdf');
    binIdx = discretize(data, edges);
    approxNLL = -sum(log(p(binIdx) + eps)); % *bw  avoid log(0)
    fprintf('bw=%.4f | approx=%.2f | true=%.2f | diff=%.4f\n', ...
        bw, approxNLL, trueNLL, approxNLL - trueNLL);
end


% plot
sigmas = [5 10 20 30 40 50];
trueNLLs = zeros(numel(sigmas), 1);
estimatedNLLs = zeros(numel(sigmas), 6);

for i=1:numel(sigmas)

    % 1) simulate data
    N = 1e4; sigma = sigmas(i);
    data = sigma * randn(N,1);
    data = mod(data + 90, 180) - 90;
    
    % 2) true NLL (analytic Gaussian)
    L = 180; K = 5;
    ks = -K:K;
    pdf_wrap = @(x) sum(normpdf(x + ks*L, 0, sigma), 2);
    trueNLL = -sum(log(pdf_wrap(data) + eps));
    % trueNLL = -sum(log(normpdf(data,0,sigma)));
    
    % save data
    trueNLLs(i) = trueNLL;

    % KDE
    [f, xi] = ksdensity(data, 'Function', 'pdf', 'Support', [-90 90]);
    likelihood = interp1(xi, f, data, 'linear'); 
    nllKDE = -log(likelihood);
    nllKDE = sum(nllKDE(:));
    fprintf("NLL KDE: %.4f True NLL: %.4f diff: %.4f\n", ...
        nllKDE, trueNLL, nllKDE - trueNLL)
    
    % 3) binned NLL for different bin widths
    % NLL artifically decreases because of binSize term in pdf
    bws = [0.02 1 2 5 10 50];
    for bw_idx = 1:numel(bws)
        bw = bws(bw_idx);
        edges = -90:bw:90; edges(end+1)=max(edges(end),90);
        p = histcounts(data, edges, 'Normalization','pdf');
        binIdx = discretize(data, edges);
        approxNLL = -sum(log(p(binIdx) + eps)); % *bw  avoid log(0)
        fprintf('bw=%.4f | approx=%.2f | true=%.2f | diff=%.4f\n', ...
            bw, approxNLL, trueNLL, approxNLL - trueNLL);

        estimatedNLLs(i, bw_idx) = approxNLL;
    end

end

figure
for i=1:numel(sigmas)
    subplot(2, 3, i)
    plot(1:numel(bws), estimatedNLLs(i, :),lineWidth=1.5)
    hold on
    yline(trueNLLs(i), lineStyle="--", lineWidth=1.5, DisplayName="true NLL")
    hold off
    ylabel("NLL")
    xlabel("binwidths")
    xticks(1:numel(bws))
    xticklabels(bws)
    legend
    title(sprintf("True data sigma %d", sigmas(i)))

    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.YTickMode = 'auto';
    ax.YTickLabelMode = 'auto';

end

