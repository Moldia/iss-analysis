%% unsupervised clustering (kmeans) of spatial data
%  transcript counts in every grid is normalized by the maximum counts
%  Xiaoyan, 2015-8-5

%% parameters
grid_size = 400;
background_image = 'input_example\IHC_mergedRGB_20%.tif';
scale = 0.2;    % image scale
image = imread(background_image);
decoded_file = 'input_example\QT_0.4_0.001_details.csv';
num_clusters = 2;

%% transcripts
tic
disp('Start processing..');
[name,Pos] = getinsitudata_f(decoded_file,2,1,0);
[name_uni,~,idx_re] = unique(name);

% remove NNNN
idx_NNNN = strcmp(name_uni,'NNNN');
if nnz(idx_NNNN)==0
    warning('Check variable name_uni.');
end
Pos = Pos(idx_re~=find(idx_NNNN),:);
idx_re = idx_re(idx_re~=find(idx_NNNN));

%% count reads in grid
grid_nrx =  ceil(max(Pos(:,1))/grid_size);
grid_nry =  ceil(max(Pos(:,2))/grid_size);
grid_nr = grid_nrx*grid_nry;

% scaling
Pos_pixel = round(correctcoordinates_f(Pos,1/grid_size));

blobinpolygon = cell(1,grid_nr);
Counted = false(length(Pos),1);

for j = 1:grid_nrx   % j along x axis, column
    fprintf('%.2f%s\n',double((j-1)/grid_nrx)*100,'% grid processed.');
    
    for i = 1:grid_nry  % i along y axis, row
        temp_in = Pos_pixel(:,1)==j & Pos_pixel(:,2)==i;
        temp_in = temp_in & ~Counted;
            
        % blobs within grid (logical)
        k = (j-1)*grid_nry+i; % counting direction: y
        blobinpolygon(k) = {temp_in};
        
        % already counted blobs
        Counted = Counted | temp_in;

    end
end
disp('100.00% grid processed.');

%% record grid positions
polygoncoord = zeros(5,2,grid_nr);
grid_size = ceil(grid_size*scale);
for j = 1:grid_nrx   % j along x axis, column
    gridxmin = (j-1)*grid_size;
    gridxmax = j*grid_size;
    
    for i = 1:grid_nry  % i along y axis, row
        gridymin = (i-1)*grid_size;
        gridymax = i*grid_size;
        
        k = (j-1)*grid_nry+i; % counting direction: y
        
        polyx = [gridxmin,gridxmin,gridxmax,gridxmax,gridxmin];
        polyy = [gridymin,gridymax,gridymax,gridymin,gridymin];
        polygoncoord(:,:,k) = [polyx',polyy'];
    end
end
clear gridxmin gridxmax gridymin gridymax polyx polyy

%% transcripts excluding NNNN
count_transcript = zeros(length(name_uni(~idx_NNNN)),grid_nr);
list = 1:length(name_uni);
% tic
for i = 1:grid_nr
    count_transcript(:,i) = (hist(idx_re(logical(blobinpolygon{i})),list(~idx_NNNN)))';
end
% toc

emptygrid = sum(count_transcript,1)==0;
count_transcript_nozero = count_transcript(:,~emptygrid);
count_transcript_norm = count_transcript_nozero./repmat(max(count_transcript_nozero,[],2),1,size(count_transcript_nozero,2));
% count_transcript_norm = count_transcript_nozero./repmat(sum(count_transcript_nozero,1),size(count_transcript_nozero,1),1);

%% figure
% figure;
% bh = bar3(count_transcript_norm);
% for i = 1:length(bh)
%     temp = get(bh(i),'ZData');
%     set(bh(i),'CData',temp);
% end
% 
% L = linkage(count_transcript_norm','average');
% set(0,'RecursionLimit',1000);
% figure; dendrogram(L)
% 
%% kmeans
disp('Start kmeans clustering with 100 replicates.');
[cidx,ctrs] = kmeans(count_transcript_norm',num_clusters,'Distance','sqeuclidean','Replicates',10);
toc

list = 1:grid_nr;

figure; 
imshow(image);
% axis image;
hold on;
col = {'red' 'green' 'blue' 'yellow' 'cyan'};
for k = 1:max(cidx)
    listtemp = list(~emptygrid);
    listtemp = listtemp(cidx==k);
    for i = 1:length(listtemp)
        temppoly = polygoncoord(:,:,listtemp(i));
        fill(temppoly(:,1),temppoly(:,2),col{k},'facealpha',0.3);
    end
end
drawnow;

figure;
bh = bar(ctrs');
set(gca,'XTick',1:length(name_uni(~idx_NNNN)),'XTIckLabel',name_uni(~idx_NNNN),...
    'XLim',[0 length(name_uni(~idx_NNNN))+1],'XTickLabelRotation',90);
for i = 1:length(bh)
    set(bh(i),'facecolor',rgb(col{i}),'edgecolor',[.2 .2 .2],'linewidth',0.1);
end
box off
title('centroid locations of clusters')
legend(gca,'show')