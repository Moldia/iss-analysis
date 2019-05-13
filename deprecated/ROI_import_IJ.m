%% import polygons from text file, plot transcripts, count and output bar plot
% example file:
% roi1
% 15864,18384,19272,19032,17808,16128,14208,13032,12720,13752,15720
% 18360,18576,20208,22440,23472,23688,23568,21888,20160,19104,18360
% 
% necrotic
% 8784,8712,9408,10560,11328,11376,9504
% 17640,16152,16128,16752,17472,18048,18096
% 
%  Xiaoyan 2015-8-26

clear; 
close all;
drawnow;

%% parameters
coordinate_file = 'input_example\1_IJ_coordinates.txt';
decoded_file = 'input_example\1_QT_0.45_1e-05_details.csv';
background_image = 'input_example\1_base2_c1_ORG_10%.tif';  % important for size
scale = 0.1;    % image scale
taglist = ID_list_TBproject_3bases; % used for plotting

name_exclude = {'mACTB'}; % exclude from plotting and bar plot
col = {'white' 'yellow' 'blue' 'tomato' 'palevioletred' 'bisque' 'pink' 'mediumspringgreen'};

show_image = 1; % for plotting
output_count_filename = 'count_ROI.csv';

%% transcripts
[name,Pos] = getinsitudata_f(decoded_file);

% unique transcripts
[name_uni,~,idx_re] = unique(name);
[p,q] = hist(idx_re,unique(idx_re));
[~,idx_sort] = sort(p,'descend');

Pos_scaled = correctcoordinates_f(Pos,scale);

%% ROI information
[ROI_number,Coord] = getIJROIinfo_f(coordinate_file);

%% count transcripts within polygons
[ROI_count,ROI_freq,ROI_proportion,ROI_area,In] = ...
    countinROI_f(ROI_number,Coord,Pos,name_uni,idx_re,p);

%% format taglist
taglist = formattaglist_f(taglist);

if show_image
    image = imread(background_image);
end

%% exclude NNNN and the specified
idx_NNNN = find(strcmp(name_uni,'NNNN'));
idx_rest = 1:length(name_uni);
if isempty(name_exclude)
    idx_rest(idx_NNNN) = [];
    idx_sort(idx_sort==idx_NNNN) = [];
else
    idx_exclude = ismember(name_uni,name_exclude);
    idx_rest(unique([idx_NNNN;find(idx_exclude)])) = [];
    idx_sort(ismember(idx_sort,[idx_NNNN;find(idx_exclude)])) = [];
end

%% bar plot
if length(col)<ROI_number
    col = repmat(col,1,ceil(ROI_number/length(col)));
end

figure;
bh = bar(ROI_freq(idx_rest,:));
set(gca,'XTick',1:length(idx_rest),'XTickLabel',name_uni(idx_rest),...
    'XLim',[0 length(idx_rest)+1],'XTickLabelRotation',90);
legend(Coord(1,:),'location','NorthEastOutside','color',[.6 .6 .6]);
ylabel('relative frequency');
box off;
for i = 1:length(bh)
    set(bh(i),'facecolor',rgb(col{i}),'edgecolor',[.2 .2 .2],'linewidth',0.1);
end
drawnow;

%% plotting
pos_sub = Pos_scaled(ismember(idx_re,idx_rest),:);
In_sub = In(ismember(idx_re,idx_rest),:);
re_idx_sub = idx_re(ismember(idx_re,idx_rest));
figure;
if show_image
	imshow(image,[]);
end
hold on;
set(gca,'YDir','reverse');
axis image;
name_legend = [];
for i = 1:length(idx_rest)
    if ~isempty(find(strcmp(taglist(:,2),name_uni(idx_rest(i)))))
        temp_pos_sub = pos_sub(logical((re_idx_sub==idx_rest(i)).*In_sub),:);
        if ~isempty(temp_pos_sub)
            name_legend = [name_legend,name_uni(idx_rest(i))];
            plot(temp_pos_sub(:,1),...
                temp_pos_sub(:,2),...
                taglist{strcmp(taglist(:,2),name_uni(idx_rest(i))),3});
        end
    end
end
legend(name_legend,'location','NorthEastOutside','color',[.6 .6 .6]);
axis off;
drawnow;

%% plotting - score
% if show_quality_score_plotting
%     score_transformed_plot = floor(score(logical(In.*ismember(idx_re,idx_rest)))/.1);
%     idx_re_plot = idx_re(logical(In.*ismember(idx_re,idx_rest)));
%     [~,~,idx_re_plot] = unique(idx_re_plot);
%     position_plot = Pos_scaled(logical(In.*ismember(idx_re,idx_rest)),:);
%     % symbols for different transcripts
%     symbol = {'o','+','*','x','s','d','^','v','>','<','p','h'};
%     if length(symbol)<length(name_uni)
%         symbol = repmat(symbol,1,ceil(length(name_uni)/length(symbol)));
%     end
%     % colors for quality score
%     col_score = parula(6);
%     plotscore_f(name_uni(idx_rest),score_transformed_plot,...
%         position_plot,...
%         idx_re_plot,...
%         show_image,background_image,symbol,col_score)
% end
% drawnow;

%% polygons
figure;
if show_image
	imshow(image,[]);
end
hold on;
set(gca,'YDir','reverse');
axis image;
for i = 1:ROI_number
    plot(Coord{2,i}*scale,Coord{3,i}*scale,'linewidth',2,'color',rgb(col{i}));
end
legend(Coord(1,:),'location','NorthEastOutside','color',[.6 .6 .6]);
axis off;
drawnow;

%% write output file
writecountfile_f(output_count_filename,1:ROI_number,Coord(1,:),name_uni,...
    ROI_count,ROI_freq,ROI_proportion,ROI_area,col)

%% clustering
% beta test, only compatible with >=R2014b
if ROI_number >= 3
    L = linkage(ROI_freq(idx_sort,:)','average');
    figure;
    sh1 = subplot('Position',[.05 .15 .05 .8]);
    dendrogram(L,'Orientation','left')
    set(sh1,'color','none',...
        'XLim',[-.001 inf],'YLim',[.5 ROI_number+.5]);
    idx_cluster = get(sh1,'YTickLabel');
    axis off
    idx_cluster = str2num(idx_cluster);
    idx_cluster = flipud(idx_cluster);
    
    sh2 = subplot('Position',[.1 .15 .8 .8]);
    bh = bar3(1:ROI_number,ROI_freq(idx_sort,idx_cluster)',1);
    
    for i = 1:length(bh)
        Zdata = get(bh(i),'ZData');
        set(bh(i),'CData',Zdata);
    end
    YLab = Coord(1,:);
    YLab = YLab(idx_cluster);
    view(2)
    set(sh2,'color','none',...
        'XLim',[.5 length(idx_rest)+.5],'YLim',[.5 ROI_number+.5],...
        'plotboxaspectratiomode','auto',...
        'YAxisLocation','Right',...
        'YTick',1:ROI_number,'YTickLabel',YLab,...
        'XTick',1:length(idx_rest),'XTickLabel',name_uni(idx_sort),...
        'XTickLabelRotation',90);
end