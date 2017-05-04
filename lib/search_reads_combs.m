function uPositives = search_reads_combs(name, pos, distlim,...
    mincomb, searchname, plotouterbox, col)
% find where any at least mincomb of different reads occur together
% Xiaoyan, 2017

uNames = unique(name);
if nargin >= 5
    [name, pos] = removereads(name, setdiff(uNames, searchname), pos);
end

% return if no reads left
if isempty(name); uPositives =[]; return; end;

[~, ~, idxName] = unique(name);

% range search
[idxNN, dNN] = rangesearch(pos, pos, distlim);
maxNN = max(cellfun(@length, idxNN));
idxNN = cellfun(@(v) [v, zeros(1, maxNN-length(v))], idxNN, 'uni', 0);
idxNN = reshape([idxNN{:}], maxNN, [])';
dNN = cellfun(@(v) [v, nan(1, maxNN-length(v))], dNN, 'uni', 0);
dNN = reshape([dNN{:}], maxNN, [])';

% step distance
positives = [];
d = distlim/10;
if nargin <= 3
    mincomb = 2;
end
while d <= distlim
    tempidx = idxNN.*(double(dNN>d-distlim/20 & dNN<=d));
    nameNN = tempidx;
    nameNN(nameNN~=0) = idxName(nameNN(nameNN~=0));
    nNN = sum(logical(nameNN), 2);
    nNNtemp = find(nNN >= mincomb);
    nameNN = nameNN(nNNtemp,:);
    positive = false(size(nameNN,1), 1);
    for i = 1:size(nameNN,1)
        uNNs = unique(nameNN(i,:));
        positive(i) = length(uNNs(uNNs~=0)) >= mincomb;
    end
    positive = nNNtemp(positive);
    
    for i = 1:length(positive)
        groups = tempidx(positive(i),:);
        groups(groups==0) = [];
        centroid = mean(pos(groups,:),1);
        positives = [positives;...
            centroid, d, pos(positive(i),:), length(unique(idxName(groups)))];
    end
    d = d + distlim/10;
end

if ~isempty(positives)
    % merge close clusters
    idxNNcluster = rangesearch(positives(:,1:2), positives(:,1:2), distlim);
    alreadyCounted = false(size(positives,1),1);
    uPositives = [];
    
    for i = 1:numel(idxNNcluster)
        idxCluster = idxNNcluster{i};
        
        % skip if the query has been processed
        if alreadyCounted(idxCluster(1)); continue; end
        
        replicates = positives(idxCluster,:);
        if size(replicates,1) > 1
            [~, sortDist] = sort(replicates(:,3));
            replicates = replicates(sortDist,:);
            nReads = replicates(:,6);
            [~, sortnReads] = sort(nReads, 'descend');
            replicates = replicates(sortnReads,:);
        end
        uPositives = [uPositives; replicates(1,:)];
        alreadyCounted(idxNNcluster{i}) = true;
    end
    
    % each query point takes only the cluster that has shortest distance
    % and most number of differerent reads
    [~, sortDist] = sort(uPositives(:,3));
    uPositives = uPositives(sortDist,:);
    nReads = uPositives(:,6);
    [~, sortnReads] = sort(nReads, 'descend');
    uPositives = uPositives(sortnReads,:);
    [~, idx] = unique(uPositives(:,4:5), 'rows');
    uPositives = uPositives(idx,:);
    
    % visualization
    if nargin <= 5
        plotouterbox = 1;
    end
    
    if nargin <= 6
        col = 'y';
    end
    
    % visualization
    boxsize = distlim*2;
    hold on;
    for i = 1:size(uPositives,1)
        plot([uPositives(i,1)-uPositives(i,3), uPositives(i,1)-uPositives(i,3), uPositives(i,1)+uPositives(i,3), uPositives(i,1)+uPositives(i,3), uPositives(i,1)-uPositives(i,3)],...
            [uPositives(i,2)+uPositives(i,3), uPositives(i,2)-uPositives(i,3), uPositives(i,2)-uPositives(i,3), uPositives(i,2)+uPositives(i,3), uPositives(i,2)+uPositives(i,3)],...
            'color', col);
        if plotouterbox
            plot([uPositives(i,4)-boxsize, uPositives(i,4)-boxsize, uPositives(i,4)+boxsize, uPositives(i,4)+boxsize, uPositives(i,4)-boxsize],...
                [uPositives(i,5)+boxsize, uPositives(i,5)-boxsize, uPositives(i,5)-boxsize, uPositives(i,5)+boxsize, uPositives(i,5)+boxsize],...
                'w', 'linewidth', 2);
        end
    end
else
    uPositives = [];
    disp('No occurence of specified reads combination within given distance.');
end

end
