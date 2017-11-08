% unsupervised clustering (kmeans) of spatial data
% transcript counts in every grid is normalized by the maximum counts
% Xiaoyan, 2017

%% input
grid_size = 400;
decoded_file = 'QT_0.35_0.0001_details_900196_1.csv';
num_clusters = 3;
output_directory = '';

background_image = '';    % can be empty if not needed for visualization
scale = 0.2;    % image scale

%% 
% transcripts
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);
[uNames, ~, iName] = unique(name);

% count reads in grid
nx =  ceil(max(pos(:,1))/grid_size);
ny =  ceil(max(pos(:,2))/grid_size);
nGrids = nx*ny;

% scaling
posGrid = round(correctcoord(pos, 1/grid_size));

% count reads in each grid
blobinpolygon = cell(1,nGrids);
Counted = false(length(pos), 1);
for j = 1:nx   % j along x axis, column
    fprintf('%.2f%s\n', double((j-1)/nx)*100, '% grid processed.');
    
    for i = 1:ny  % i along y axis, row
        ingrid = posGrid(:,1)==j & posGrid(:,2)==i;
        ingrid = ingrid & ~Counted;
            
        % blobs within grid (logical)
        k = (j-1)*ny+i; % counting direction: y
        blobinpolygon(k) = {ingrid};
        
        % already counted blobs
        Counted = Counted | ingrid;

    end
end
disp('100.00% grid processed.');

% record grid positions
cx = (meshgrid(1:nx, 1:ny)-1)*grid_size + grid_size/2;
cy = (meshgrid(1:ny, 1:nx)'-1)*grid_size + grid_size/2;
cPolygons = sortrows([cx(:), cy(:)]);

coordPolygons = zeros(5, 2, nGrids);
szGrid = ceil(grid_size*scale);
for j = 1:nx   % j along x axis, column
    gridxmin = (j-1)*szGrid;
    gridxmax = j*szGrid;
    
    for i = 1:ny  % i along y axis, row
        gridymin = (i-1)*szGrid;
        gridymax = i*szGrid;
        
        k = (j-1)*ny+i; % counting direction: y
        
        polyx = [gridxmin,gridxmin,gridxmax,gridxmax,gridxmin];
        polyy = [gridymin,gridymax,gridymax,gridymin,gridymin];
        coordPolygons(:,:,k) = [polyx',polyy'];
    end
end
clear gridxmin gridxmax gridymin gridymax polyx polyy

% gene counts
cGenes = zeros(length(uNames), nGrids);
for i = 1:nGrids
    cGenes(:,i) = (hist(iName(logical(blobinpolygon{i})), 1:numel(uNames)))';
end

emptygrid = sum(cGenes,1)==0;
nonEmptyGrids = find(~emptygrid);
cGenes_nonzero = cGenes(:,~emptygrid);

% write files (all data)
cNames = {'grid_num', uNames{:}, 'center_x', 'center_y'};

fid = fopen(fullfile(output_directory, 'GridClustering_GeneCount_all.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids; cGenes_nonzero; cPolygons(nonEmptyGrids,:)']);
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

% select genes for clustering
cbValues = checkboxes(uNames);
idx = cellfun(@(v) find(strcmp(v, uNames)), cbValues(:,1));
isSelected = false(numel(idx), 1);
isSelected(idx) = cell2mat(cbValues(:,2));

cGenes_selected = cGenes_nonzero(isSelected,:);
emptygrid = sum(cGenes_selected,1)==0;
nonEmptyGrids = nonEmptyGrids(~emptygrid);
cGenes_selected = cGenes_selected(:,~emptygrid);

% normalized by max
cGenes_maxnorm = bsxfun(@rdivide, cGenes_selected, max(cGenes_selected,[],2));

% normalization by sum
cGenes_sumnorm = bsxfun(@rdivide, cGenes_selected, sum(cGenes_selected,1));


% kmeans
disp('Starting kmeans clustering with 100 replicates..');
[iCluster, centroid] = kmeans(cGenes_maxnorm', num_clusters,...
    'Distance', 'sqeuclidean', 'Replicates', 100);

% write files (only selected)
cNames = {'grid_num', uNames{isSelected}, 'center_x', 'center_y'};

fid = fopen(fullfile(output_directory, 'GridClustering_GeneCount.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids; cGenes_selected; cPolygons(nonEmptyGrids,:)']);
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

fid = fopen(fullfile(output_directory, 'GridClustering_GeneCount_SumNorm.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids; cGenes_sumnorm; cPolygons(nonEmptyGrids,:)']);
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

cNames = {'grid_num', uNames{isSelected}, 'center_x', 'center_y', 'cluseter_id'};
fid = fopen(fullfile(output_directory, 'GridClustering_GeneCount_MaxNorm.csv'), 'w');
fprintf(fid, lineformat('%s', numel(cNames)), cNames{:});
towrite = num2cell([nonEmptyGrids; cGenes_maxnorm; cPolygons(nonEmptyGrids,:)'; iCluster']);
fprintf(fid, lineformat('%d', numel(cNames)), towrite{:});
fclose(fid);

% visualization
figure; 
try
    image = imread(background_image);
    imshow(image);
catch
    axis image
    set(gca, 'YDir', 'reverse');
end
col = {'red' 'green' 'blue' 'yellow' 'cyan'};
hold on;
for k = 1:max(iCluster)
    grids = nonEmptyGrids(iCluster==k);
    for i = 1:length(grids)
        temppoly = coordPolygons(:,:,grids(i));
        fill(temppoly(:,1), temppoly(:,2), col{k}, 'facealpha', 0.3);
    end
end
drawnow;

% heatmap of normalized counts and barplot of cluster centroid position
figure;
ax1 = subplot(121);
[~, idxSort] = sort(iCluster);
[~, idxFirst] = unique(iCluster(idxSort));
imagesc(cGenes_maxnorm(:,idxSort));
hold on;
plot(repmat(idxFirst', 2, 1), repmat([0; nnz(isSelected)+1], 1, numel(idxFirst)),...
    'r', 'linewidth', 2);
xlabel('bin');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', uNames(isSelected), 'fontsize', 5);
title('normailzed bin count data')
colorbar

ax2 = subplot(122);
bh = barh(centroid');
set(gca, 'ytick', 1:nnz(isSelected), 'yticklabel', uNames(isSelected),...
    'ylim',[0 nnz(isSelected)+1], 'fontsize', 5, 'ydir', 'reverse');
for i = 1:length(bh)
    set(bh(i), 'facecolor', rgb(col{i}), 'edgecolor', [.2 .2 .2], 'linewidth', .1);
end
box off
xlabel('normalized gene count')
title('centroid location of clusters')
legend(catstrnum('Cluster ', 1:max(iCluster)))

linkaxes([ax1, ax2], 'y');
