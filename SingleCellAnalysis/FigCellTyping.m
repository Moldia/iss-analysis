% demo figure for in situ cell typing
% Xiaoyan, 2016-10-16
% data: 160408_160222UCL_1-hippo1


%% cell clustering
cellblob = importdata('..\CP_160511_NeuN\Stitched\Nuclei_blobs.csv');
genes = cellblob.textdata(1,2:end-2);
cellblob = cellblob.data; 
pos = cellblob(:,end-1:end);
cellid = cellblob(:,1);
Data = cellblob(:,2:end-2);

LData = Data./repmat(max(Data,[],1),length(Data),1);    % relative scale

% hierarchical
L = linkage(LData);
D = dendrogram(L,0);
idx_cell = get(gca,'XTickLabel');
idx_cell = str2num(idx_cell);
idx_cell = flipud(idx_cell);
C = cluster(L,'maxclust',10);

L2 = linkage(LData');
D2 = dendrogram(L2,0);
idx_gene = get(gca,'XTickLabel');
idx_gene = str2num(idx_gene);
idx_gene = flipud(idx_gene);

figure;
bh = bar3(Data(idx_cell,idx_gene),1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
    set(bh(i),'edgecolor','none');
end
view(2)
axis normal
set(gca,'xtick',1:length(idx_gene),...
    'xticklabel',genes(idx_gene),'xticklabelrotation',90,...
    'ylim',[.5 length(idx_cell)+.5]);

% % tSNE
% T = tsne(Data);
% figure,plot(T(:,1),T(:,2),'.');

%% top-level classification visualization
xrange = [3655.55582129418,4067.23461002292];
yrange = [11485.6464992877,11958.2736473318];
genes = {
    'Plp1'
    'Rorb'
    'Calb1'
    'Nrn1'
    'Reln'
    'Penk'
    'Slc6a1'
    'Lhx6'
    'Ndnf'
    'Cxcl14'
    'Cnr1'
    'Vip'
    'Htr3a'
    'Chodl'
    'Sst'
    'Pvalb'
};
colors = {    
    'blue'
    'red'
    'red'
    'red'
    'green'
    'green'
    'green'
    'green'
    'green'
    'green'
    'green'
    'green'
    'green'
    'green'
    'green'
    'green'
};

load('..\CP_160511_NeuN\Stitched\RenumberedCells.mat');
LabelImg = plottogether_f('..\CP_160511_NeuN\Stitched\Nuclei_blobs.csv',...
    '..\base1_c1_ORG.tif',Irenumbered,genes,colors);
set(gca,'xlim',xrange,'ylim',yrange);

figure; imshow(LabelImg); set(gca,'xlim',xrange,'ylim',yrange);

%% with original data
uiopen('..\CP_160424\Plotted.fig',1);
hold on;
h = imshow(LabelImg);
set(h,'alphadata',.2);
set(gca,'xlim',xrange,'ylim',yrange);

%% interneuron subtype visualization
genes = {'Pvalb','Vip','Cabl2','Htr3a','Reln','Sst','Cxcl14'};
colors = {'purple','orange','orange','orange','cyan','cyan','red'};

LabelImg = plottogether_f('..\CP_160511_NeuN\Stitched\Nuclei_blobs.csv',...
    '..\base1_c1_ORG.tif',Irenumbered,genes,colors);
set(gca,'xlim',xrange,'ylim',yrange);

