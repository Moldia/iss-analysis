% proximity analysis
% rewritten from proximity_analysis_main.R and proximity_analysis_source.R
% original R scripts by Olle
% rewritten by Xiaoyan, 2014-9-28
% Xiaoyan 2014-11-18


clear; 
close all;

%% parameters to set
decoding_file = 'Decoding\QT_0.3_0.01_details.csv';
pair_distance = 100;
name_subset_exclude = {'NNNN' 'ACTB'};  % transcripts to exclude in NN analysis
number_neighbors = 50;
rep = 50;   % number of replicates
subpair = 10;   % the number of pairs to show in the submatrix (recommended not to change it)

%% extract data
% import
[name,pos] = getinsitudata_f(decoding_file);

% data subset
name = name(~ismember(name,name_subset_exclude));
pos = pos(~ismember(name,name_subset_exclude),:);

if number_neighbors+1>size(pos,1)
    error(['Dataset is too small to find ' num2str(number_neighbors) ' of NNs.']);
end

%% calculations
[Score_neighbor,Score_probability,Score_dev,name_uni,name_p,num_uni] = ...
    neighborknn_f(name,pos,number_neighbors,pair_distance,rep);

%% sub-heatmap with most significant pairs
[a,b] = sort(Score_dev(:),'descend');
sub_idx = b(1:subpair*2);
sub_idx = sub_idx(mod(1:subpair*2,2)~=0);
n = ceil(sub_idx/num_uni);
m = mod(sub_idx,num_uni);
m(m==0)=num_uni;

Score_dev_submatrix = Score_dev(m,n);

figure;
bh = bar3(1:subpair,Score_dev_submatrix);
for i = 1:length(bh)
     zdata = get(bh(i),'Zdata');
     set(bh(i),'Cdata',zdata)
end
view(2); colorbar;
set(gca,'yDir','rev','xAxisLocation','top',...
    'xtick',1:subpair,'ytick',1:subpair,...
    'xticklabel',name_uni(n),'yticklabel',name_uni(m));
axis([0 subpair+1 0 subpair+1]);

set(gcf,'name','deviation heatmap sub-matrix(sorted)');
xticklabel_rotate([],90)

%% table to show the abundances of transcripts shown in the submatrix
f=figure;
set(f,'units','normalized','position',[0.05 0.2 0.2 0.4],...
    'name','count table(submatrix)');
h = uitable(f,'data',[name_uni(m),num2cell(name_p(m))',name_uni(n),num2cell(name_p(n))'],...
    'ColumnName',{'name' 'count' 'name' 'count'},...
    'ColumnFormat',{'char' 'numeric' 'char' 'numeric'},...
    'RowName',1:subpair);
set(h,'units','normalized','position',[0 0 1 1]);

%% hierarchichal clustering
% Dist = ones(size(Score_probability))./abs(Score_neighbor);
% Dist(Dist==Inf) = 10;
% Z = linkage(Dist);
% figure;dendrogram(Z)
% temp = get(gca,'XTickLabel');
% temp = str2num(temp);
% set(gca,'XTickLabel',name_uni(temp));
% set(gca,'XTIckLabelRotation',90)
% 
% figure;
% bh = bar3(1:length(name_uni),Score_dev(temp,temp));
% for i = 1:length(bh)
%     zdata = get(bh(i),'Zdata');
%     set(bh(i),'Cdata',zdata)
% end
% view(2); colorbar;
% tick = 1:length(name_uni);
% tick = tick(mod(tick,5)==0);
% name_sort = name_uni(temp);
% tick_sort = name_sort(mod(1:num_uni,5)==0);
% 
% set(gca,'yDir','rev','xAxisLocation','top',...
%     'xtick',tick,'ytick',tick,...
%     'xticklabel',tick_sort,'yticklabel',tick_sort);
% axis([0 length(name_uni)+1 0 num_uni+1]);
% set(gcf,'name',figure_name);