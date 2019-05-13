
[name, pos, quality] = getinsitudata('CP_1\beforeQT_details.csv');

% find neighbors
[nnIdx, nnDist] = rangesearch(pos, pos, 5);
nNeighbor = cellfun(@length, nnIdx) - 1;
pNNNN = zeros(max(nNeighbor)+1,1);
for i = 0:max(nNeighbor)
    subreads = nNeighbor==i;
    pNNNN(i+1) = nnz(strcmp('NNNN', name(subreads)))/nnz(subreads);
end

ax1 = subplot(211);
distributionPlot(quality, 'groups', nNeighbor)
set(gca, 'xticklabel', 0:max(nNeighbor))
xlabel('number of neighbors within 5 px radius');
ylabel('quality score');

subplot(212);
[ax2, lax1, lax2] = plotyy(0:max(nNeighbor), pNNNN,...
    0:max(nNeighbor), log2(hist(nNeighbor, 0:max(nNeighbor))));
lax1.LineWidth = 1.5;
lax2.LineWidth = 1.5;
box off
set(ax2(1), 'xtick', 0:max(nNeighbor), 'xlim', ax1.XLim);
set(ax2(2), 'xtick', 0:max(nNeighbor), 'xlim', ax1.XLim);
ylabel(ax2(1), 'proportion of NNNN')
ylabel(ax2(2), 'log2 of read count')
    
% distributionPlot(strcmp(name,'NNNN'),...
%     'groups', nNeighbor, 'histOpt', 0, 'showMM', 0);
% set(gca, 'ylim', [-.5 1.5], 'ytick', [0 1], 'yticklabel', {'expected', 'NNNN'});

