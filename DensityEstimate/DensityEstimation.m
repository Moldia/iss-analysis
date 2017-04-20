% kernel density estimation of a specific gene
% Xiaoyan, 2017

clear;
close all;
drawnow;

%% parameters
decoded_file = 'D:\Salamander project\Ch030\CP_170418\Decoding\QT_0.55_0.0001_details_ROI.csv';
image = 'D:\Salamander project\Ch030\Preprocessing\Stitched\Ch030_170310_Pw_5dpa_SBL2_CX_c1_stitched.tif';   % important for size
scale = 1;      % image scale
name_density = '60S';
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
