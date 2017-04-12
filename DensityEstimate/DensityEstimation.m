% kernel density estimation of a specific gene
% Xiaoyan, 2017

clear;
close all;
drawnow;

%% parameters
decoded_file = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\spots_ROI1.csv';
image = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\Ab_c1_ROI1.tif';    % important for size
scale = 1;      % image scale
name_density = 'Pvalb';
bandwid = 50;   % in original scale

%%
% all transcripts
[name, pos] = getinsitudata(decoded_file);
pos = correctcoord(pos, .2);

% unique transcripts
[uniName, ~, idxName] = unique(name);
[p, q] = hist(idxName, 1:length(uniName));

% image size
img = imfinfo(image);
imsize = [img.Height, img.Width];

% density estimation plot
idx = find(strcmp(uniName, name_density));
if isempty(idx)
    error('No specified transcript detected in the input file');
end
pos_density = pos(idxName==idx, :);

if size(pos_density,1)>2
    [~, density] = kde2d_modified(pos_density, 2^10, [0 0],...
        floor([imsize(2)/scale/5 imsize(1)/scale/5]), floor([bandwid/5 bandwid/5]));
else
    density = zeros(2^10);
end
density = imresize(density, [imsize(1)/5 imsize(2)/5]);

figure;
imshow(density,[]);
colormap(gca, parula);
title(name_density);