function [NamesToShow, ColorsToShow] = pie_plot_size_white(o)
% plot pie chart for each cell, showing probability of it belonging to all
% classes

nC = size(o.CellYX,1);
nK = size(o.pCellClass,2);

% find classes to collapse
CollapseMe = zeros(nK,1);
Colors = zeros(nK,3);
DisplayName = o.ClassNames;
for i=1:size(o.ClassCollapse,1)
    ClassList = o.ClassCollapse{i,1};
    for j=1:length(ClassList)
        MyClasses = strmatch(ClassList{j}, o.ClassNames);
        if length(MyClasses)==0;  continue; end
        CollapseMe(MyClasses)=i;
        Colors(MyClasses,:) = repmat(o.ClassCollapse{i,3},length(MyClasses),1);
        DisplayName(MyClasses) = o.ClassCollapse(i,2);
    end
end

nColorWheel = sum(CollapseMe==0);

Colors0 = hsv(ceil(nColorWheel*1.2));
Colors(~CollapseMe,:) = Colors0(1:nColorWheel,:); % last is zero

figure(43908765)
% figure
clf; 
set(gcf, 'Color', 'w');
set(gca, 'color', 'w');
hold on

% load(o.CellMapFile, 'RelCellRadius');

for c=1:nC
   
    pMy = o.pCellClass(c,:);
    
%     % merge three subtypes of PC.CA1
%     pMy(1) = sum(pMy(1:3));
%     pMy(2:3) = 0;
    
    % sum up probabilities of merged classes
    for j = 1:max(CollapseMe)
        collapseThese = find(CollapseMe==j);
        pMy(collapseThese(1)) = sum(pMy(CollapseMe==j));
        pMy(collapseThese(2:end)) = 0;
    end
    
    
    
    WorthShowing = find(pMy>o.MinPieProb);
    if ~isempty(WorthShowing)

        h = pie(pMy(WorthShowing), repmat({''}, 1, sum(WorthShowing>0)));

        for i=1:length(h)/2
            hno = (i*2-1);
    %         Index = find(strcmp(h(i*2).String, NickNames), 1);
            set(h(hno), 'FaceColor', Colors(WorthShowing(i),:));
            
%             % all same size
%             set(h(hno), 'Xdata', get(h(hno), 'Xdata')*o.PieSize + o.CellYX(c,2));
%             set(h(hno), 'Ydata', get(h(hno), 'Ydata')*o.PieSize + o.CellYX(c,1));

%             % size based on relative cell radius
%             set(h(hno), 'Xdata', get(h(hno), 'Xdata')*o.PieSize*RelCellRadius(c) + o.CellYX(c,2));
%             set(h(hno), 'Ydata', get(h(hno), 'Ydata')*o.PieSize*RelCellRadius(c) + o.CellYX(c,1));            

            % size based on number of reads
            set(h(hno), 'Xdata', get(h(hno), 'Xdata')*o.PieSize*sum(o.pSpotCell(:,c)) + o.CellYX(c,2));
            set(h(hno), 'Ydata', get(h(hno), 'Ydata')*o.PieSize*sum(o.pSpotCell(:,c)) + o.CellYX(c,1));            
            
%             set(h(hno), 'EdgeAlpha', 0);
            set(h(hno), 'EdgeAlpha', 1, 'LineWidth', .1, 'EdgeColor', [.5 .5 .5]);
        end
    end
    
    if mod(c,2000)==0
        drawnow
    end
end

yMax = max(o.CellYX(:,1));
xMax = max(o.CellYX(:,2));
yMin = min(o.CellYX(:,1));
xMin = min(o.CellYX(:,2));

ClassShown = find(any(o.pCellClass>o.MinPieProb,1));
ClassDisplayNameShown = DisplayName(ClassShown);
[uDisplayNames, idx] = unique(ClassDisplayNameShown, 'stable');
nShown = length(uDisplayNames);

% % sort based on color
% [uCols, ~, legend_order] = unique(Colors(ClassShown((idx)),[2,1,3]), 'rows');
% [~, legend_order] = sort(legend_order);
legend_order = [1:numel(idx)]';

% move PC, non neuron and uncalled to end
bottom = [find(strcmp(uDisplayNames, 'Uncalled'))
    find(strcmp(uDisplayNames, 'Non neuron'))
    find(strcmp(uDisplayNames, 'PC other'))
    find(strcmp(uDisplayNames, 'PC CA1'))
    find(strcmp(uDisplayNames, 'Calb2.Cryab'))];
top = find(strcmp(uDisplayNames, 'Sst.Cryab'));
legend_order = [bottom; setdiff(legend_order, [bottom; top], 'stable');...
    top];

NamesToShow = DisplayName(ClassShown(idx(legend_order)));
ColorsToShow = Colors(ClassShown(idx(legend_order)),:);
for k=1:nShown
    h = text(xMax*1.1 - xMin*.2, yMin + k*(yMax-yMin)/nShown, NamesToShow{k}, 'fontsize', 8);
    set(h, 'color',ColorsToShow(k,:));
end

