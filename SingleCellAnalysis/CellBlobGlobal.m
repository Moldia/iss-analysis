%% relate blobs to cells
%  single cell analysis
%  Xiaoyan, 2015-7-6


clear;

%% renumber cells
tiletoremove = [1,2,3,11,12,17,28,29,31,44,59,106,121,134,136,137,138,148:165];

cells = csvread('Cells.csv',1);
cellstart = 0;
tileuni = unique(cells(:,3));
for i = 1:length(tileuni)
    if ismember(tileuni(i),tiletoremove)
        cells(cells(:,3)==tileuni(i),8) = 0;    % renumber identified cells in empty tiles to 0
    else
        cells(cells(:,3)==tileuni(i),8) = cells(cells(:,3)==tileuni(i),2)+cellstart;
        cellstart = cellstart + max(cells(cells(:,3)==tileuni(i),2));
    end
end

%% cell global position
tile_start_pos = getposition('E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\ZENout\Tiled_150702_CellBlob.csv');

cellswrite = [cells(:,8),cells(:,3),cells(:,5:6)];
tilemax = max(cellswrite(:,2));
cellpos = zeros(size(cellswrite,1),2);

for t = 1:tilemax
    if ~isempty(cellswrite(cellswrite(:,2)==t,:))
        cellpos(cellswrite(:,2)==t,:) = bsxfun(@plus,cellswrite(cellswrite(:,2)==t,3:4),tile_start_pos(tile_start_pos(:,1)==t,2:3));
    end
end

cellswrite = [cellswrite,cellpos];

fid = fopen('cells_renumbered.csv','w');    
fprintf(fid,'CellID,TileID,Tile_X_pos,Tile_Y_pos,Global_X_pos,Global_Y_pos\n');
fprintf(fid,'%d,%d,%d,%d,%d,%d\n',cellswrite');
fclose(fid);

%% merge information
cellblobs = csvread('blobs.csv',1);

relation = cellblobs(:,[2,3,6]);

blobstart = 0;

tileuni = unique(relation(:,2));

tic
for i = 1:length(tileuni)
    blobsintile = relation(:,2)==tileuni(i);
    blobtile = relation(blobsintile,:);

    blobtile(:,1) = blobtile(:,1) + blobstart;   % blob ID
    parentcelltile = blobtile(:,3);   % parent cells

    [cell_uni,~,idx_re] = unique(parentcelltile);
    celltile = cells(cells(:,3)==tileuni(i),:);
    
    for j = 1:length(cell_uni)
        if cell_uni(j)~=0
            cellid = celltile(celltile(:,7)==cell_uni(j),8);
            blobtile(idx_re==j,3) = cellid;
        end
    end
    blobstart = blobstart + length(relation(relation(:,2)==tileuni(i)));
    relation(blobsintile,:) = blobtile;
end
toc

%% summarize after QT, noNNNN
load('E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\CP_150624\Decoding\QT_0.4_0.001.mat');
relation = [relation,global_x_pos,global_y_pos];
relation_allbt = relation(blob_allbt,:);
relation_allqt = relation(blob_allqt,:);

fid = fopen('signal_cell.csv','w');
fprintf(fid,'blob_ID,tile_ID,parent_cell,global_x_pos,global_y_pos\n');
fprintf(fid,'%d,%d,%d,%d,%d\n',relation');
fclose(fid);

fid = fopen('signal_cell_allbt.csv','w');
fprintf(fid,'blob_ID,tile_ID,parent_cell,global_x_pos,global_y_pos\n');
fprintf(fid,'%d,%d,%d,%d,%d\n',relation_allbt');
fclose(fid);

fid = fopen('signal_cell_allqt.csv','w');
fprintf(fid,'blob_ID,tile_ID,parent_cell,global_x_pos,global_y_pos\n');
fprintf(fid,'%d,%d,%d,%d,%d\n',relation_allqt');
fclose(fid);

% noNNNN subsets
relation_allqt_noNNNN = relation_allqt;
relation_allqt_noNNNN(strcmp(name_tag_allqt,'NNNN'),:) = [];
name_tag_allqt_noNNNN = name_tag_allqt;
name_tag_allqt_noNNNN(strcmp(name_tag_allqt,'NNNN'),:) = [];
allqt_noNNNN = allqt;
allqt_noNNNN(strcmp(name_tag_allqt,'NNNN'),:) = [];

% transform data
[~,~,tag_order] = unique(expected_list);
[transformed_data,c] = hist3([relation_allqt_noNNNN(:,3),allqt_noNNNN],{unique(relation_allqt_noNNNN(:,3)),unique(allqt_noNNNN)});
transformed_data = [(unique(expected_list))';transformed_data];
transformed_data = transformed_data(:,tag_order);
sum(transformed_data(2:end,:),1)

% merge same names
transformed = [];
[unitag,first_idx,re_idx] = unique(exp_tags,'stable');
for i = 1:length(unitag)
    transformed = [transformed,sum(transformed_data(2:end,re_idx==i),2)];
end
cell_blob = [c{1,1}',transformed];
cell_transformed_data = [{'cellID'},unitag;num2cell(cell_blob)];
sum(transformed,1)

fid = fopen('cell_signal.csv','w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',cell_transformed_data{1,:});
fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',cell_blob');
fclose(fid);
