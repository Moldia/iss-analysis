function plotscore_f(name_uni,score_transformed,position,idx_re,...
    show_image,background_image,symbol,col_matrix)
% colors for quality score
% symbols for different transcripts
% Xiaoyan, 2014-12-16

figure;
subplot(1,10,1:9);
if show_image
    image = imread(background_image);
    imshow(image,[]);
else
    axis image;
    set(gca,'YDir','reverse');
end
axis off;
hold on;

temp_score = find(ismember(score_transformed,6:10));
temp = position(temp_score,:);
temp_idx_re = idx_re(temp_score);
for i = 1:length(name_uni)
    position_plot = temp(temp_idx_re==i,:);
    if ~isempty(position_plot)
        plot(position_plot(:,1),position_plot(:,2),...
            'LineStyle','none','Marker',symbol{i},...
            'MarkerFaceColor',col_matrix(6,:),'MarkerEdgeColor',col_matrix(6,:));
    end
end

for j = 5:-1:0
    temp_score = find(score_transformed==j);
    temp = position(temp_score,:);
    temp_idx_re = idx_re(temp_score);
    
    if ~isempty(temp_score)
        for i = 1:length(name_uni)
            position_plot = temp(temp_idx_re==i,:);
            if ~isempty(position_plot)
                plot(position_plot(:,1),position_plot(:,2),...
                    'LineStyle','none','Marker',symbol{i},...
                    'MarkerFaceColor',col_matrix(j+1,:),'MarkerEdgeColor',col_matrix(j+1,:));
            end
        end
    end
end
          

subplot(1,10,10);
hold on;
for j = 0:5
    plot(1,j,'o','MarkerFaceColor',col_matrix(j+1,:),'MarkerEdgeColor',col_matrix(j+1,:)');
end
for j = 6:10
    plot(1,j,'o','MarkerFaceColor',col_matrix(6,:),'MarkerEdgeColor',col_matrix(6,:));
end
set(gca,'YTick',0:10,'YTickLabel',0:.1:1,'XTick',[]);
title({'Color scheme' 'for quality scores'});