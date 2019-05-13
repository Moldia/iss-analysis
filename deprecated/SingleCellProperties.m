
[name, pos, quality] = getinsitudata('CP_1\beforeQT_details.csv');

% find neighbors
[nnIdx, nnDist] = rangesearch(pos, pos, 10);
nNeighbor = cellfun(@length, nnIdx) - 1;
pNNNN = zeros(max(nNeighbor)+1,1);
for i = 0:max(nNeighbor)
    subreads = nNeighbor==i;
    pNNNN(i+1) = nnz(strcmp('NNNN', name(subreads)))/nnz(subreads);
end

subplot(211);
distributionPlot(quality, 'groups', nNeighbor)
set(gca, 'xticklabel', 0:max(nNeighbor))
xlabel('number of neighbors within 10 px radius');
ylabel('quality score');

set(gca, 'xtick', 0:max(nNeighbor))
ax2 = axes('Position',ax1.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
plot(0:max(nNeighbor), pNNNN, 'parent', ax2);
linkaxes([ax1, ax2], 'xy')
    
    
distributionPlot(strcmp(name,'NNNN'),...
    'groups', nNeighbor, 'histOpt', 0, 'showMM', 0);
set(gca, 'ylim', [-.5 1.5], 'ytick', [0 1], 'yticklabel', {'expected', 'NNNN'});

