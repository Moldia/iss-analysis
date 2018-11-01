% plot semi-transparent circles for selected subgroup of transcritps
% Xiaoyan, 2017

close all;

%% modify here
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
output_folder = 'E:\PROOOJECTS\test_dataset\TransparentPlot'; 
image = 'E:\PROOOJECTS\test_dataset\860502_1_align.png';
scale = .2;
names = {'EpCAM', 'CDH1'};

show_image = 1;
radius = 200;       % original scale
transparency = .2;      % between 0 and 1. 1 is not not transparent 

%% do not modify

% load
[name, pos] = getinsitudata(decoded_file);
pos = correctcoord(pos, scale);

% unique transcripts
[uNames, ~, iNames] = unique(name);

% circles
radius = radius*scale;
alph = linspace(0, 2*pi, 15);
circx = radius*cos(alph);
circy = radius*sin(alph);

% plot
figure;
if show_image
    img = imread(image);
    imshow(img);
else
    plotonblank;
end

hold on;

try
    mymap = [1 0 0; 0 1 0; 0 0 1; 1 1 1];
    col = mymap(1:length(names),:);
catch
    col = lines(length(names));
end
figure_title = [];
for i = 1:length(names)
    idx = find(strcmp(uNames, names{i}));
    if isempty(idx)
        warning(['Could not find ' names{i} ' in the input file.'])
    else
        x = pos(iNames==idx,1);
        y = pos(iNames==idx,2);
        plotx = (repmat(x,1,15)+repmat(circx,length(x),1))';
        ploty = (repmat(y,1,15)+repmat(circy,length(y),1))';
        patch(plotx, ploty, col(i,:), 'facealpha', transparency, 'linestyle', 'none');
        figure_title = [figure_title, names(i)];
    end
end
title(strjoin(figure_title, '-'));
