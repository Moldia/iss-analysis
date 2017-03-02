%% Overlapping circles to visualize neighbors
%  for Eph data
%  transcripts from the different section

%%
format compact
warning('off','all');

%% parameters
image = 'E:\PROOOJECTS\9_Ephrin\Sample_snapshot\8128_3_5_s1c1.jpg';
scale = 1;

% all transcripts defined below will be shown as RED
decoded_file_list_1 = {...    % a list of files to read
    'E:\PROOOJECTS\9_Ephrin\Image analysis\8128_A2_EGFR\8128_3_1\Decoded_details.csv',...
    'E:\PROOOJECTS\9_Ephrin\Image analysis\8128_A2_EGFR\8128_5_1\Decoded_details.csv'};
channel_order_list_1 = {...   % original channel, four channels for each file
    'EGFR mut' 'A2 wt' 'EGFR wt' 'A2 mut';...
    'EGFR mut' 'A2 wt' 'EGFR wt' 'A2 mut'};
name_list_1 = {...
    'EGFR mut';...
    'EGFR mut'}; % transcripts to plot, one for each section

% all transcripts defined below will be shown as GREEN
decoded_file_list_2 = {...    % a list of files to read
    'E:\PROOOJECTS\9_Ephrin\Image analysis\8128_A3_A5\8128_3_2\Decoded_details.csv',...
    'E:\PROOOJECTS\9_Ephrin\Image analysis\8128_A3_A5\8128_5_2\Decoded_details.csv'};
channel_order_list_2 = {...   % original channel, four channels for each file
    'A3 mut' 'A5 mut' 'A5 wt' 'A3 wt';...
    'A3 mut' 'A5 mut' 'A5 wt' 'A3 wt'};
name_list_2 = {...
    'A5 mut';...
    'A5 mut'}; % transcripts to plot, one for each section

radius = 100;
transparent_level = 0.3;    % between 0 and 1
figure_title = 'EGFR - A5';
show_image = 1;

%% extract data - different sections
Data1 = struct;
for s = 1:length(decoded_file_list_1)
    data = csvread(decoded_file_list_1{s},1);
    field = ['file' num2str(s)];
    Data1.(field).position = data(:,1:2);
    Data1.(field).uniname = channel_order_list_1(s,:);
    Data1.(field).uniname_re = data(:,4);
end
Data2 = struct;
for s = 1:length(decoded_file_list_2)
    data = csvread(decoded_file_list_2{s},1);
    field = ['file' num2str(s)];
    Data2.(field).position = data(:,1:2);
    Data2.(field).uniname = channel_order_list_2(s,:);
    Data2.(field).uniname_re = data(:,4);
end

%% reform to pair info
pos_1 = [];
for s = 1:length(decoded_file_list_1)
    field = ['file' num2str(s)];
    position = Data1.(field).position;
    uniname = Data1.(field).uniname;
    uniname_re = Data1.(field).uniname_re;
    
    idx = find(strcmp(uniname,name_list_1{s}));
    if isempty(idx)
        warning(['Could not find ' name_list_1{s} ' in the input file.'])
    else
        x = position(uniname_re==idx,1);
        y = position(uniname_re==idx,2);
        pos_1 = [pos_1;x,y];
    end
end
pos_2 = [];
for s = 1:length(decoded_file_list_2)
    field = ['file' num2str(s)];
    position = Data2.(field).position;
    uniname = Data2.(field).uniname;
    uniname_re = Data2.(field).uniname_re;
    
    idx = find(strcmp(uniname,name_list_2{s}));
    if isempty(idx)
        warning(['Could not find ' name_list_2{s} ' in the input file.'])
    else
        x = position(uniname_re==idx,1);
        y = position(uniname_re==idx,2);
        pos_2 = [pos_2;x,y];
    end
end

%% overlapping
[plotx_1,ploty_1] = dot2poly_f(pos_1(:,1),pos_1(:,2),radius,15);
[plotx_2,ploty_2] = dot2poly_f(pos_2(:,1),pos_2(:,2),radius,15);

% find intersections
[I,~] = rangesearch(pos_2,pos_1,2*radius);   % first find neighbors within a fixed radius to save computational time

NN_num = zeros(length(I),1);
for i = 1:length(I)
    NN_num(i) = size(I{i},2);
end

idx_query = find(NN_num);
plotx_inter = {}; ploty_inter = {};

for i = 1:length(idx_query)
    idx_pool = I{idx_query(i)};
    for j = 1:length(idx_pool)
        [x,y] = polybool('intersection',...
            plotx_1(:,idx_query(i)),ploty_1(:,idx_query(i)),...
            plotx_2(:,idx_pool(j)),ploty_2(:,idx_pool(j)));
        
        if isempty(x)
        else
            plotx_inter = [plotx_inter; {x}];
            ploty_inter = [ploty_inter; {y}];
        end
    end
end

%% plot
figure; hold on;
if show_image
    I = imread(image);
    imshow(I,[]);
end
axis image;
axis off;
set(gca,'YDir','reverse');

col = {'r' 'g'};

patch(plotx_1*scale,ploty_1*scale,col{1},'facealpha',transparent_level,'linestyle','none');
patch(plotx_2*scale,ploty_2*scale,col{2},'facealpha',transparent_level*1.2,'linestyle','none');

for i = 1:length(plotx_inter)
    patch(plotx_inter{i}*scale,ploty_inter{i}*scale,'y','facealpha',max(1,2*transparent_level),'linestyle','none');
end
title(figure_title);


