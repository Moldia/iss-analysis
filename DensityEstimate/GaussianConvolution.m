%% density plot (guassian smoothing)
%  Xiaoyan, 2015-4-7

clear;
% close all;
drawnow;

%% parameters
decoded_file = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\spots_ROI1.csv';
image = 'K:\161230_161220KI_3-1\Kenneth\hippocampi\hippocampus1\Ab_c1_ROI1.tif';    % important for size
name_density = 'Pvalb';
bandwid = 50;

%% transcripts
[name,Pos] = getinsitudata(decoded_file,1,1);

% unique transcripts
[name_uni,~,idx_re] = unique(name);
[p,q] = hist(idx_re,unique(idx_re));

%% image size
imgin = imfinfo(image);
Isize = [imgin.Height,imgin.Width];

%% gaussian smoothing
idx_gaussian = find(strcmp(name_uni,name_density));
if isempty(idx_gaussian)
    error('No specified transcript detected in the input file');
end
pos_gaussian = Pos(idx_re==idx_gaussian,1:2);

temp = floor(pos_gaussian/5);
temp(temp==0) = 1;
Itemp = accumarray(fliplr(temp),1,floor(Isize/5));
fh = fspecial('gaussian',bandwid*2,bandwid/5);
Itemp = imfilter(Itemp,fh);

figure;
imshow(Itemp/max(fh(:)),[]);
colormap(hot);
colorbar
title(name_density);