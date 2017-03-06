% overlay density estimates of two genes
% Xiaoyan, 2017

clear;
close all;
drawnow;

%% parameters
decoded_file = 'Decoding\beforeQT_details.csv';
image = 'AlignedDAPI_b1.png';    % important for size
scale = .2;
name_density = {'CC10', 'Drd3'};       % first in red, second in cyan
bandwid = 50;  % in original scale
use_absolute_count = 0;

%% 
% all transcripts
[name,pos] = getinsitudata(decoded_file);
pos = correctcoord(pos,.2);

% unique transcripts
[uniName, ~, idxName] = unique(name);
[p, q] = hist(idxName, 1:length(uniName));

% image size
img = imfinfo(image);
imsize = [img.Height, img.Width];

% density estimation plot
idx = cellfun(@(v) strcmp(v, uniName), name_density, 'uni', 0);
try
    idx = cellfun(@find, idx);
catch
    error('At least one of the genes specified not found in input file.');
end

I = zeros(floor(imsize(1)/5), floor(imsize(2)/5), 3, 'double');
layer = [1 0 0; 0 1 1];
for i = 1:2
    pos_density = pos(idxName==idx(i),:);
    if use_absolute_count
        temp = floor(pos_density/5);
        temp(temp==0) = 1;
        Itemp = accumarray(fliplr(temp), 1, floor(imsize/5));
        fh = fspecial('gaussian', bandwid*2, bandwid/5);
        smoothed = imfilter(Itemp, fh);
        smoothed = smoothed/max(fh(:));     % normalize to peak height of gaussian filter
    else
        if size(pos_density,1)>2
            [bandwidth, smoothed, X, Y] = kde2d_modified(pos_density, 2^10, [0 0],...
                floor([imsize(2)/scale/5 imsize(1)/scale/5]), floor([bandwid/5 bandwid/5]));
        else
            smoothed = zeros(2^10);
        end
        smoothed = imresize(smoothed, floor([imsize(1)/5 imsize(2)/5]));
        smoothed = (smoothed - min(smoothed(:)))/max(smoothed(:));
    end
    
    I = I + cat(3, smoothed*layer(i,1), smoothed*layer(i,2), smoothed*layer(i,3));
end

I = I/max(I(:))*2;
figure;
imshow(I);
if use_absolute_count
    title([name_density{1} ' - ' name_density{2} ' gaussian smoothed' ]);
else
    title([name_density{1} ' - ' name_density{2} ' densities' ]);
end

