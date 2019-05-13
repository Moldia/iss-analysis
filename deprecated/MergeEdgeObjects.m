%% merge objects identified multiple times due to edge effect
%  Xiaoyan, 2015-7-2
% deprecated


close all;
clear;

%%
tilesize_x = 2000;
tilesize_y = 2000;
tilenum = [11,15];

%% shifted image
load('E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\CP_150702_Cell_Shifted\DefaultOut.mat');

cellpos_y = handles.Measurements.Nuclei.Location_Center_Y;
cellpos_x = handles.Measurements.Nuclei.Location_Center_X;

cellpos = [];
imagenum = [];
tile = [];
objnum = [];
for i = 1:length(cellpos_y)
    if ~isempty(cellpos_y{i})
        imagenum = [imagenum;repmat(i,length(cellpos_y{i}),1)];
        tile = [tile;repmat(double(handles.Measurements.Image.Metadata_position{i}),length(cellpos_y{i}),1)];
        objnum = [objnum;(1:length(cellpos_y{i}))'];
        cellpos = [cellpos; [cellpos_x{i},cellpos_y{i}]];
    end
end

cells = [imagenum,objnum,tile,objnum,cellpos,objnum,objnum];
clearvars -except cells relation_allqt_noNNNN pos_allqt tilenum tilesize_x tilesize_y


cellstart = 0;
tileuni = unique(cells(:,3));
for i = 1:length(tileuni)
    cells(cells(:,3)==tileuni(i),8) = cells(cells(:,3)==tileuni(i),2)+cellstart;
    cellstart = cellstart + max(cells(cells(:,3)==tileuni(i),2));
end

% global position
tile_start_pos = getposition('E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\ZENout\Tiled_150702_Cell_Shifted.csv');

cellswrite = [cells(:,8),cells(:,3),cells(:,5:6)];
tilemax = max(cellswrite(:,2));
cellpos = zeros(size(cellswrite,1),2);

for t = 1:tilemax
    if ~isempty(cellswrite(cellswrite(:,2)==t,:)) && ismember(t,tile_start_pos(:,1))
        cellpos(cellswrite(cellswrite(:,2)==t,1),:) = bsxfun(@plus,cellswrite(cellswrite(:,2)==t,3:4),tile_start_pos(tile_start_pos(:,1)==t,2:3));
    end
end

cellswrite = [cellswrite,cellpos];

cellpos_new = cellpos+1;


%% original image
load('DefaultOut.mat');

cellpos_y = handles.Measurements.Nuclei.Location_Center_Y;
cellpos_x = handles.Measurements.Nuclei.Location_Center_X;

cellpos = [];
imagenum = [];
tile = [];
objnum = [];
for i = 1:length(cellpos_y)
    if ~isempty(cellpos_y{i})
        imagenum = [imagenum;repmat(i,length(cellpos_y{i}),1)];
        tile = [tile;repmat(double(handles.Measurements.Image.Metadata_position{i}),length(cellpos_y{i}),1)];
        objnum = [objnum;(1:length(cellpos_y{i}))'];
        cellpos = [cellpos; [cellpos_x{i},cellpos_y{i}]];
    end
end

cells = [imagenum,objnum,tile,objnum,cellpos,objnum,objnum];
clearvars -except cells cellpos_new relation_allqt_noNNNN pos_allqt tilenum tilesize_x tilesize_y


cellstart = 0;
tileuni = unique(cells(:,3));
for i = 1:length(tileuni)
    cells(cells(:,3)==tileuni(i),8) = cells(cells(:,3)==tileuni(i),2)+cellstart;
    cellstart = cellstart + max(cells(cells(:,3)==tileuni(i),2));
end

% global position
tile_start_pos = getposition('E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\ZENout\Tiled_150702_CellBlob.csv');

cellswrite = [cells(:,8),cells(:,3),cells(:,5:6)];
tilemax = max(cellswrite(:,2));
cellpos = zeros(size(cellswrite,1),2);

for t = 1:tilemax
    if ~isempty(cellswrite(cellswrite(:,2)==t,:)) && ismember(t,tile_start_pos(:,1))
        cellpos(cellswrite(cellswrite(:,2)==t,1),:) = bsxfun(@plus,cellswrite(cellswrite(:,2)==t,3:4),tile_start_pos(tile_start_pos(:,1)==t,2:3));
    end
end

cellswrite = [cellswrite,cellpos];
cellpos = cellpos+1;

%% create plot line image and distance transform
D_threshold = 1;
Isize = imfinfo('E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\150616_4028_1_2_b1_mip_20%_c1.jpg');
scale = 0.2;
Isize = [Isize.Height,Isize.Width];
Isize = ceil(Isize/scale);
plotline = zeros(ceil(Isize/5));

for i = 1:tilenum(1)-1
    plotline(tilesize_y/5*i,:)=1;
end
for i = 1:tilenum(2)-1
    plotline(:,tilesize_x/5*i)=1;
end
plotline_dist = bwdist(plotline);
% plotline_dist = max(plotline_dist(:))-plotline_dist;

plotline_buffer = plotline_dist;
plotline_buffer(plotline_buffer==0) = 1;
plotline_buffer(plotline_buffer>D_threshold) = 0;
plotline_buffer(plotline_buffer~=0) = 1;

%% plot
cellpos_scaled = cellpos/5;
cellpos_new_scaled = cellpos_new/5;
cellpos_new_scaled(cellpos_new_scaled(:,1)<0 | cellpos_new_scaled(:,2)<0,:) = [];

hold on;
plot(cellpos_scaled(:,1),cellpos_scaled(:,2),'.');
plot(cellpos_new_scaled(:,1),cellpos_new_scaled(:,2),'.');
set(gca,'YDir','reverse');
axis image
legend

for i = 1:tilenum(1)-1
    plot([0 tilesize_x*tilenum(2)]/5,[tilesize_y*i+1,tilesize_y*i+1]/5,'b');
end
for i = 1:tilenum(2)-1
    plot([tilesize_x*i+1,tilesize_x*i+1]/5,[0 tilesize_y*tilenum(1)]/5,'b');
end

ih = imshow(plotline_buffer);
set(ih,'alphadata',0.5);

%% find nuclei/cells in the buffering zone
cellpos_scaled_linear = (ceil(cellpos_scaled(:,1))-1)*Isize(1)/5+ceil(cellpos_scaled(:,2));
buffer_linear = reshape(plotline_buffer,[],1);
cellpos_scaled_buffer = cellpos_scaled(buffer_linear(cellpos_scaled_linear)==1,:);
hold on;
plot(cellpos_scaled_buffer(:,1),cellpos_scaled_buffer(:,2),'o');

%% distance threshold for the new positions
plotline_buffer = plotline_dist;
plotline_buffer(plotline_buffer==0) = 1;
plotline_buffer(plotline_buffer>1) = 0;
plotline_buffer(plotline_buffer~=0) = 1;
buffer_linear = reshape(plotline_buffer,[],1);

cellpos_new_scaled_linear = (ceil(cellpos_new_scaled(:,1))-1)*Isize(1)/5+ceil(cellpos_new_scaled(:,2));
cellpos_new_scaled_buffer = cellpos_new_scaled(buffer_linear(cellpos_new_scaled_linear)==1,:);
plot(cellpos_new_scaled_buffer(:,1),cellpos_new_scaled_buffer(:,2),'o');

%% convolve
convolution = zeros(ceil(Isize/5));
convolution = reshape(convolution,[],1);
convolution(cellpos_scaled_linear(buffer_linear(cellpos_scaled_linear)==1)) = 1;
temp = zeros(ceil(Isize/5));
temp = reshape(temp,[],1);
temp(cellpos_new_scaled_linear(buffer_linear(cellpos_new_scaled_linear)==1)) = 1;
convolution = convolution+temp;
convolution = reshape(convolution,Isize(1)/5,[]);

st = strel('disk',2);
convolution = imfilter(convolution,double(st.getnhood));
figure, imshow(convolution,[]);
colormap jet

area_recover = convolution >= 3;
% imshow(area_recover,[]);
area_recover = imfilter(area_recover,double(st.getnhood));
imshow(area_recover,[]);

%% nuclei in the recover area
conn_area_recover = bwconncomp(area_recover);
L = double(labelmatrix(conn_area_recover));
L = L(:);

celltorecover = ismember(cellpos_scaled_linear,find(area_recover));
cellspot_linear = cellpos_scaled_linear(celltorecover);

hold on;
plot(cellpos_scaled(celltorecover,1),cellpos_scaled(celltorecover,2),'.');
merged_cell_id = L(cellspot_linear);

celltorecover = find(celltorecover);

%% visualize
relation_allqt = csvread('signal_cell_allqt.csv',1);

figure;
imshow('E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\150616_4028_1_2_b1_mip_20%_c1.jpg');
hold on;
sym = {'*','o','^','h','p'};
col = hot(max(merged_cell_id));

for i = 1:max(merged_cell_id)    
    celltomerge = celltorecover(merged_cell_id==i);
    for j = 1:length(celltomerge)
        blobs = find(ismember(relation_allqt(:,3),celltomerge(j)));
        if ~isempty(blobs)
            plot(relation_allqt(blobs,4)/5,relation_allqt(blobs,5)/5,'color',col(i,:),'marker',sym{j});
            relation_allqt(blobs,3) = celltomerge(1);
        end
    end
end

fid = fopen('signal_cell_allqt_edgecorrected.csv','w');
fprintf(fid,'blob_ID,tile_ID,parent_cell,global_x_pos,global_y_pos\n');
fprintf(fid,'%d,%d,%d,%d,%d\n',relation_allqt');
fclose(fid);

