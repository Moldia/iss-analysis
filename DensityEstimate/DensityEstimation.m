% kernel density estimation of a specific gene
% Xiaoyan, 2017

clear;
close all;
drawnow;

%% parameters
decoded_file = 'C:\Worky\Temp\QT_0.4_1e-20_details.csv';
image = 'C:\Worky\Temp\TissueA_Fluo_staining.png';   % important for size
scale = 1;      % image scale
name_density = 'HER2';
bandwid = 50;   % in original scale

%%
% all transcripts
[name, pos] = getinsitudata(decoded_file);

% density estimation plot
density = gene_kde(name, pos, name_density, bandwid, image, scale);

figure;
imshow(density, []);
colormap(gca, parula);
title(name_density);
