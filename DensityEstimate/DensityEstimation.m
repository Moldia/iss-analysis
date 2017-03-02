%% kernel density estimation of a specific gene
%   Xiaoyan 2015-8-10

clear;
close all;
drawnow;

%% parameters
decoded_file = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\spots_ROI1.csv';
image = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\Ab_c1_ROI1.tif';    % important for size
scale = 1;
name_density = 'Pvalb';
bandwid = 50;  % in original scale

%% transcripts
[name,Pos] = getinsitudata(decoded_file,1,1);
Pos = correctcoord(Pos,.2);

% unique transcripts
[name_uni,~,idx_re] = unique(name);
[p,q] = hist(idx_re,unique(idx_re));

%% image size
imgin = imfinfo(image);
Isize = [imgin.Height,imgin.Width];

%% density estimation plot
idx_density = find(strcmp(name_uni,name_density));
if isempty(idx_density)
    error('No specified transcript detected in the input file');
end
pos_density = Pos(idx_re==idx_density,1:2);

if size(pos_density,1)>2
    [bandwidth,density,X,Y]=kde2d_modified(pos_density,2^10,[0 0],floor([Isize(2)/scale/5 Isize(1)/scale/5]),floor([bandwid/5 bandwid/5]));
else
    density = zeros(2^10);
end
density = imresize(density,[Isize(1)/5 Isize(2)/5]);

figure;
imshow(density,[]);
colormap(parula);
title(name_density);