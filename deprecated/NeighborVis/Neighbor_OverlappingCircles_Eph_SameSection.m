%% Overlapping circles to visualize neighbors
%  for Eph data
%  transcripts from the same section
%  red: the first transcript in the figure title
%  green: the second transcript in the figure title
%  yellow: overlapping

%%
format compact
warning('off','all');

%% parameters
decoded_file = 'E:\PROOOJECTS\9_Ephrin\Image analysis\8128_A2_EGFR\8128_3_1\Decoded_details.csv';
image = 'E:\PROOOJECTS\9_Ephrin\Sample_snapshot\8128_3_5_s1c1.jpg';
scale = 1;
channel_order = {'EGFR mut' 'A2 wt'	'EGFR wt' 'A2 mut'};    % original channel order

pair = {'EGFR wt' 'EGFR mut'};
radius = 100;
transparent_level = 0.3;    % between 0 and 1

show_image = 1;

%% transcripts
data = csvread(decoded_file,1);
pos = data(:,1:2);

% unique transcripts
name_uni = channel_order;
idx_re = data(:,4);

% pair gene index
name_p1 = find(strcmp(name_uni,pair{1}));
name_p2 = find(strcmp(name_uni,pair{2}));
if isempty(name_p1) || isempty(name_p2)
    error('At least one of the genes specified does not have any positional information.');
end

%% overlapping
[x1,y1,x2,y2,interx,intery,pool_id] = ...
    pairintersection_f(name_p1,name_p2,idx_re,pos,radius);

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
if pool_id == 2
    col = fliplr(col);
end

patch(x1*scale,y1*scale,col{1},'facealpha',transparent_level,'linestyle','none');
patch(x2*scale,y2*scale,col{2},'facealpha',transparent_level*1.2,'linestyle','none');

for i = 1:length(interx)
    patch(interx{i}*scale,intery{i}*scale,'y','facealpha',max(1,2*transparent_level),'linestyle','none');
end
title([pair{1} ' - ' pair{2} ' (r=' num2str(radius) ')']);


