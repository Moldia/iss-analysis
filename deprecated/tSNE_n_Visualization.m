% some clustering on binned or single-cell data
% input format: rows-cells/bins, columns-genes
% ONLY works in version >=2018a!!!!!
% Xiaoyan, 2018

clear;
close all;

%% modify here
hexbin_file = 'E:\PROOOJECTS\test_dataset\HexbinClustering_GeneCount_300.csv';
hexbin_size = 300; % if single-cell data, hexbin_size = 0;
remove_genes = {'ACTB'};
output_directory = 'E:\PROOOJECTS\test_dataset\tsne\';

%% do not modify

% load
data = importdata(hexbin_file);

genes = data.colheaders(:,2:end-2);
pos = data.data(:,end-1:end);
counts = data.data(:,2:end-2);
counts(:,ismember(genes, remove_genes)) = [];
genes = setdiff(genes, remove_genes);

% PCA
[coeff, score, latent] = pca(counts);
% visualize first two components
figure, biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', genes);
title('top two principle components');

% tSNE in MATLAB
% ONLY >=R2018a
seeds = 1e-4*randn(size(counts,1), 3);
Y = tsne(counts, 'NumDimensions', 3, 'NumPCAComponents', 50, 'Perplexity', 30,...
    'Standardize', 1, 'LearnRate', 1000, 'Verbose', 1, 'InitialY', seeds); 
figure, plot(Y(:,1),Y(:,2), '.');
title('tSNE dim reduction to three, shown first two');

% visualize tSNE in RGB
if ~hexbin_size;	hexbin_size = 5;    end
figure; hold on;
scale = 1;
Yrgb = rgbscale(Y);
for i = 1:size(Yrgb,1)
    [gridR, gridL, xypos] = hexbin_coord2grid(pos(i,:), hexbin_size);
    [vy, vx] = dot2poly(pos(i,2)*scale,pos(i,1)*scale,...
        hexbin_size*scale, 6);
    patch(vx, vy, Yrgb(i,:));
    %     plot(pos(idx,1), pos(idx,2), 'o',...
    %         'color', Yrgb(i,:),...
    %         'linewidth', 3);
end

% write
mkdir(output_directory);
csvwrite(fullfile(output_directory, 'tSNE_3D.csv'), Y);
csvwrite(fullfile(output_directory, 'tSNE_initial.csv'), seeds);



