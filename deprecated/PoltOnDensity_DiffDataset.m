%% plotting on top of the density estimate
%  Xiaoyan, 2015-1-10

clear;
close all;
drawnow;

%% parameters
decoded_file = 'input_example\QT_0.3_0.007_details.csv';
taglist = ID_Jessica_2;

image = 'input_example\Nuclei870276.tif';    % original full size image, important for size
decoded_file_density = 'input_example\QT_0.3_0.007_details.csv';
name_density = 'COL3A1';
bandwid = 600;

%% transcripts
[name,Pos] = getinsitudata_f(decoded_file);
position = Pos/5;

% unique transcripts
[name_uni,~,idx_re] = unique(name);
[p,q] = hist(idx_re,unique(idx_re));

%% image size
imgin = imfinfo(image);
Isize = [imgin.Height,imgin.Width];

%% prepare density estimation plot
[name2,Pos] = getinsitudata_f(decoded_file_density);
% unique transcripts
[name_uni2,~,idx_re2] = unique(name2);

idx_density = find(strcmp(name_uni2,name_density));
if isempty(idx_density)
    error('No specified transcript detected in the input file');
end
pos_density = Pos(idx_re2==idx_density,1:2);
if size(pos_density,1)>2
    [bandwidth,density,X,Y]=kde2d_X(pos_density,2^10,[0 0],[Isize(2) Isize(1)],[bandwid bandwid]);
else
    density = zeros(2^10);
end
density = imresize(density,[Isize(1)/5 Isize(2)/5]);

figure;
imshow(density,[]);
colormap(parula);
hold on;

%% plotting
taglist = formattaglist_f(taglist);
name_plot = taglist(:,2);
leg = name_plot;
for i = 1:length(name_plot)
    idx = find(strcmp(name_uni,name_plot{i}));
    if isempty(idx)
        warning(['Could not find ' name_plot{i} ' in the input file.']);
        leg(i) = [];
    else
        x = position(idx_re==idx,1);
        y = position(idx_re==idx,2);
        plot(x,y,taglist{i,3});
    end
end
title(['density - ' name_density]);
legend(leg)