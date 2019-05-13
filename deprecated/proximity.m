% proximity analysis
% rewritten from proximity_analysis_main.R and proximity_analysis_source.R
% original R scripts by Olle
% rewritten by Xiaoyan, 140725
% check proximity_2 for an updated version


clear all;


%=====parameters=====
pair_distance = 10;
name_subset_exclude = {'NNNN' 'ACTB'};
number_neighbors = 15;
rep = 100;  % the number of replicates


%---import data---
data_raw = importdata('proximity_analysis_sample_data_A549_celline.csv',',',1);
name = data_raw.textdata(2:end,2);

% take data subset
name_subset = name(~ismember(name,name_subset_exclude));
data_subset = data_raw.data(~ismember(name,name_subset_exclude),:);
pos_subset = data_subset(:,1:2);

% unique transcripts
[name_uni,idx_first,idx_re] = unique(name_subset);
[name_p,name_q] = hist(idx_re,1:length(name_uni));
[name_count_sorted,idx_sort] = sort(name_p,'descend');

%---kd nearest neighbor search---
score_neighbor = zeros(length(name_uni));
blobnum = 1:size(pos_subset(:,1));

for m = 1:length(name_uni);
    % query and pool position lists
    query_set = pos_subset(strcmp(name_subset,name_uni(m)),:);
    pool = pos_subset(~strcmp(name_subset,name_uni(m)),:);
    idx_pool = blobnum(~strcmp(name_subset,name_uni(m)));
    [idx_neighbor,d] = knnsearch(pool,query_set,...
        'k',min(number_neighbors,size(pool,1)),...
        'NSmethod','kdtree',...
        'Distance','euclidean');
    field = ['name' num2str(m)];
    Dist.(field) = d;   % distance of number_neighbors nearest neighbors
    idx_neighbor = idx_pool(idx_neighbor);
    Neighbor.(field) = idx_neighbor;    % the original indeces of neighbors
    
    % distance thresholding
    d(d>pair_distance) = 0;   
    DistThres.(field) = d;  % distance of neighbors after thresholding
    
    idx_neighbor(d==0) = 0;
    idx_neighbor = reshape(idx_neighbor',[],1);
    NeighborThres.(field) = idx_neighbor(idx_neighbor~=0);  % indices of neighbors closer than the distance threshold
    
    neighborname = idx_re(idx_neighbor(idx_neighbor~=0));
    NeighborNameThres.(field) = neighborname;   % names of the nearest neighbors closer than the distance threshold
    if length(unique(neighborname))==1
        q = neighborname(1);
        p = length(neighborname);     
    else
        [p,q] = hist(neighborname,unique(neighborname));
    end
    score_neighbor(m,q) = p;    % occurance of neighbors
end



%---generate random datasets---
score_neighbor_random = zeros(length(name_uni),length(name_uni),rep);

for r = 1:rep
    disp(r);
    % bootstrap re-sampling
    randm = randperm(size(data_subset,1));
    name_subset_random = name_subset(randm,:);
    blobnum = randm;
    
    % same kd-tree nearest neighbor searching
    for m = 1:length(name_uni)
        query_set = pos_subset(strcmp(name_subset_random,name_uni(m)),:);
        pool = pos_subset(~strcmp(name_subset_random,name_uni(m)),:);
        idx_pool = blobnum(~strcmp(name_subset_random,name_uni(m)));

        [idx_neighbor,d] = knnsearch(pool,query_set,...
            'k',min(number_neighbors,size(pool,1)),...
            'NSmethod','kdtree',...
            'Distance','euclidean');
        
        idx_neighbor = idx_pool(idx_neighbor);
        d(d>pair_distance) = 0;
        
        idx_neighbor(d==0) = 0;
        idx_neighbor = reshape(idx_neighbor',[],1);
        
        neighborname = idx_re(idx_neighbor(idx_neighbor~=0));
            
        if length(unique(neighborname))==1
            q = neighborname(1);
            p = length(neighborname);
        else
            [p,q] = hist(neighborname,unique(neighborname));
        end
        score_neighbor_random(m,q,r) = p;
        if ismember(m,q)
            error(num2str(m));
        end
    end
    
end


%---probability and deviation---
score_probability = zeros(length(name_uni));
score_dev= zeros(length(name_uni));

for i = 1:length(name_uni)
    for j = 1:length(name_uni)
        if i == j
            % skip if querry = target
        else
            [f,x] = ecdf(reshape(score_neighbor_random(i,j,:),[],1));
            if score_neighbor(i,j) >= max(x) && score_neighbor(i,j)
                score_probability(i,j) = 1;
                score_dev(i,j) = (score_neighbor(i,j)-mean(reshape(score_neighbor_random(i,j,:),[],1)))/std(reshape(score_neighbor_random(i,j,:),[],1));
            elseif score_neighbor(i,j) < min(x) || score_neighbor(i,j)==0
            else
                [m,n] = min(abs(x-score_neighbor(i,j)));
                score_probability(i,j) = f(n);
                score_dev(i,j) = (score_neighbor(i,j)-mean(reshape(score_neighbor_random(i,j,:),[],1)))/std(reshape(score_neighbor_random(i,j,:),[],1));
            end  
        end
    end
end
score_count_sorted = score_probability(fliplr(idx_sort),idx_sort);
score_dev_sorted = score_dev(fliplr(idx_sort),idx_sort);

%---images---
tick = 1:length(name_uni);
tick = tick(mod(tick,5)==0);
figure;
% bh = bar3(1:length(name_uni),score_count_sorted);  % sorted matrix
bh = bar3(1:length(name_uni),score_probability);
for i = 1:length(bh)
     zdata = get(bh(i),'Zdata');
     set(bh(i),'Cdata',zdata)
end
view(2); colorbar;
% set(gca,...
%     'xtick',tick,'ytick',tick,...
%     'xticklabel',name_count_sorted(mod(1:length(name_count_sorted)-1,5) == 0),...
%     'yticklabel',fliplr(name_count_sorted(mod(1:length(name_count_sorted)-1,5) == 0)));
set(gca,...
    'xtick',tick,'ytick',tick,...
    'xticklabel',tick,...
    'yticklabel',tick);



figure;
% bh = bar3(1:length(name_uni),score_dev_sorted);
bh = bar3(1:length(name_uni),score_dev);
for i = 1:length(bh)
     zdata = get(bh(i),'Zdata');
     set(bh(i),'Cdata',zdata)
end
view(2); colorbar;
% set(gca,...
%     'xtick',tick,'ytick',tick,...
%     'xticklabel',name_count_sorted(mod(1:length(name_count_sorted)-1,5) == 0),...
%     'yticklabel',fliplr(name_count_sorted(mod(1:length(name_count_sorted)-1,5) == 0)));
set(gca,...
    'xtick',tick,'ytick',tick,...
    'xticklabel',tick,...
    'yticklabel',tick);