% use bin count data to run k-means clusetering
% Xiaoyan, 2017

%% input
bin_count_file = 'GridClustering_GeneCount_MaxNorm.csv';  % requires header
num_clusters = 3;
output_directory = '';


%%
% import data
tableCount = readtable(bin_count_file, 'ReadVariableNames', 1);

% original column names
cNames = tableCount.Properties.VariableNames;

% make sure there are no multiple entries of the same gene
assert(numel(cNames)== numel(unique(cNames)),...
    'Column names are not unique!')

% create checkboxes and get selected values
cbValues = checkboxes(cNames);
idx = cellfun(@(v) find(strcmp(v, cNames)), cbValues(:,1));
isSelected = false(numel(idx), 1);
isSelected(idx) = cell2mat(cbValues(:,2));
cGenes = table2array(tableCount(:,isSelected));

% k-means clustering
disp('Starting kmeans clustering with 100 replicates..');
[iCluster, centroid] = kmeans(cGenes, num_clusters,...
    'Distance', 'sqeuclidean', 'Replicates', 100);

% output
fid = fopen(fullfile(output_directory, 'kmeans.csv'), 'w');
% assume first column is some sort of identifier
rNames = table2cell(tableCount(:,1));
header = [{'bin_id'}, cNames(isSelected), {'cluster_id'}];
fprintf(fid, lineformat('%s', numel(header)), header{:});

if isnumeric(rNames{1})
    fmt = lineformat('%d', numel(header));
else
    fmt = ['%s', lineformat('%d', numel(header)-1)];
end
towrite = [table2cell(tableCount(:,1)),...
    num2cell(cGenes), num2cell(iCluster)]';
fprintf(fid, fmt, towrite{:});
fclose(fid);

% heatmap of normalized counts and barplot of cluster centroid position
figure;
ax1 = subplot(121);
[~, idxSort] = sort(iCluster);
[~, idxFirst] = unique(iCluster(idxSort));
imagesc(cGenes(idxSort,:)');
hold on;
plot(repmat(idxFirst', 2, 1), repmat([0; nnz(isSelected)+1], 1, numel(idxFirst)),...
    'r', 'linewidth', 2);
xlabel('bin');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', cNames(isSelected), 'fontsize', 4);
title('normailzed bin count data')
colorbar

ax2 = subplot(122);
bh = barh(centroid');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', cNames(isSelected),...
    'ylim',[0 nnz(isSelected)+1], 'fontsize', 4, 'ydir', 'reverse');
box off
xlabel('normalized gene count')
title('centroid location of clusters')
legend(catstrnum('Cluster ', 1:max(iCluster)))

linkaxes([ax1, ax2], 'y');
