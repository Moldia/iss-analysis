function [Score_neighbor,Score_probability,Score_dev,name_uni,name_p,num_uni] = ...
    neighborknn(name,pos,number_neighbors,pair_distance,rep)

% find possibly co-localizing transcipts
% pair-wise
% original script in R by Olle
% Xiaoyan 2014-11-18

%% unique transcripts
[name_uni,~,idx_re] = unique(name);
num_uni = length(name_uni);

[name_p,~] = hist(idx_re,1:num_uni);
[name_count_sorted,idx_sort] = sort(name_p,'descend');  % sort based on abundances

%% kd-tree nearest neighbor serach
blobnum = 1:size(pos(:,1));

% kd-tree nearest neighbor search
[idx_NS,Dist] = knnsearch(pos,pos,...
    'k',number_neighbors+1,...
    'NSmethod','kdtree',...
    'Distance','euclidean');
name_NS = idx_re(idx_NS); % the gene id of neighbors

% distance thresholding
Dist(Dist>pair_distance) = 0;
name_NS(Dist==0) = 0;

% warning
if ~isempty(find(Dist(:,end),1))
    warning(['There are more potential neighbors than '...
        num2str(number_neighbors) ' within the distance of '...
        num2str(pair_distance) ' pixels.']);
end

Score_neighbor = neighborscore(name_NS,idx_re);

%% boot-strapping
Score_neighbor_random = zeros(num_uni,num_uni,rep);

for r = 1:rep
    if mod(r-1,10)==0 && r~=1
        disp([num2str(r-1) ' replicates have been finished.']);
    end
    % randomize transcripts
    random_order = randperm(length(blobnum));
    random_idx_re = idx_re(random_order);
    
    % use previously saved NN search result
    randon_name_NS = random_idx_re(idx_NS);
    randon_name_NS(Dist==0) = 0;
    
    Score_neighbor_random(:,:,r) = neighborscore(randon_name_NS,random_idx_re);
    
end

%% scores and sorting
[Score_probability,Score_dev] = bootscore;
Score_probability_sorted = Score_probability(idx_sort,idx_sort);
Score_dev_sorted = Score_dev(idx_sort,idx_sort);

%% plot - unsorted
tick = 1:length(name_uni);
tick = tick(mod(tick,2)==0);

heatmap(Score_probability,tick,tick,'probability heatmap (unsorted)');
heatmap(Score_dev,tick,tick,'deviation heatmap (unsorted)');

% look up table
f=figure;
set(f,'units','normalized','position',[0.05 0.2 0.14 0.4],...
    'name','look-up table(unsorted)');
h = uitable(f,'data',[name_uni(:) num2cell(name_p')],...
    'ColumnName',{'name' 'count'},...
    'ColumnFormat',{'char' 'numeric'},...
    'RowName',1:num_uni);
set(h,'units','normalized','position',[0 0 1 1]);

%% plot - sorted
tick = 1:length(name_uni);
tick = tick(mod(tick,5)==0);
tick_sort = name_count_sorted(mod(1:num_uni,5)==0);


heatmap(Score_probability_sorted,tick,tick_sort,'probability heatmap (sorted)');
heatmap(Score_dev_sorted,tick,tick_sort,'deviation heatmap (sorted)');

% look-up table
f=figure;
set(f,'units','normalized','position',[0.05 0.2 0.14 0.4],...
    'name','look-up table(sorted)');
h = uitable(f,'data',name_uni(idx_sort),...
    'ColumnName',{'name'},...
    'ColumnFormat',{'char'},...
    'RowName',name_count_sorted);
set(h,'units','normalized','position',[0 0 1 1]);


%%
    function Score = neighborscore(name_NS,idx)
        Score = zeros(num_uni);
        
        % group into tanscript pairs
        for m = 1:num_uni
            neighbor_temp = name_NS(idx==m,:);   % gene m's neighbors
            temp = neighbor_temp(neighbor_temp~=0); % take non-zeros
            if isempty(temp)
            else
                [p,~] = hist(temp,1:num_uni);
                Score(m,:) = p;
            end
        end
    end


    function [Score_probability,Score_dev] = bootscore
        Score_probability = zeros(num_uni);
        Score_dev = zeros(num_uni);
        
        for i = 1:num_uni
            for j = 1:num_uni
                if i == j
                    % skip if querry = target
                else
                    [f,x] = ecdf(reshape(Score_neighbor_random(i,j,:),[],1));
                    if Score_neighbor(i,j) >= max(x) && Score_neighbor(i,j)
                        Score_probability(i,j) = 1;
                        Score_dev(i,j) = (Score_neighbor(i,j)-mean(reshape(Score_neighbor_random(i,j,:),[],1)))/std(reshape(Score_neighbor_random(i,j,:),[],1));
                    elseif Score_neighbor(i,j) < min(x) || Score_neighbor(i,j)==0
                    else
                        [~,n] = min(abs(x-Score_neighbor(i,j)));
                        Score_probability(i,j) = f(n);
                        Score_dev(i,j) = (Score_neighbor(i,j)-mean(reshape(Score_neighbor_random(i,j,:),[],1)))/std(reshape(Score_neighbor_random(i,j,:),[],1));
                    end
                end
            end
        end
    end


    function heatmap(plotdata,tick,tick_sort,figure_name)
        figure;
        bh = bar3(1:num_uni,plotdata);
        for i = 1:length(bh)
            zdata = get(bh(i),'Zdata');
            set(bh(i),'Cdata',zdata)
        end
        view(2); colorbar;
        
        set(gca,'yDir','rev','xAxisLocation','top',...
            'xtick',tick,'ytick',tick,...
            'xticklabel',tick_sort,'yticklabel',tick_sort);
        axis([0 length(name_uni)+1 0 num_uni+1]);
        set(gcf,'name',figure_name);
    end

end