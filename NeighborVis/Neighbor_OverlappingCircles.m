%% Overlapping circles to visualize neighbors
%  pair-wise
%  Xiaoyan 2015-1-20

%%
format compact
warning('off','all');

%% parameters
decoding_file = 'input_example\DECODE_0.4_0.007_beforeQT_details_newform.csv';
pair = {'COL3A1' 'HER2'};
radius = 300;
transparent_level = 0.15;    % between 0 and 1
image = 'input_example\RGB_adjusted3.jpg';
show_image = 1;

%% transcripts
[name,pos] = getinsitudata_f(decoding_file);
 
% unique transcripts
[name_uni,~,idx_re] = unique(name);

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

patch(x1,y1,col{1},'facealpha',transparent_level,'linestyle','none');
patch(x2,y2,col{2},'facealpha',transparent_level*1.2,'linestyle','none');

for i = 1:length(interx)
    patch(interx{i},intery{i},'y','facealpha',max(1,2*transparent_level),'linestyle','none');
end
title([pair{1} ' - ' pair{2} ' (r=' num2str(radius) ')']);


