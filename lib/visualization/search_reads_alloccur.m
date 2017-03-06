function positives = search_reads_alloccur(name, pos, distlim)
% find where all reads occur together
% Xiaoyan, 2017


% use the least abundant one as query
[uniName, ~, idxName] = unique(name);
countName = hist(idxName, 1:length(uniName));
[~, idxQ] = sort(countName, 'ascend');
idxQ = idxQ(1);
posquery = pos(idxName==idxQ,:);

% range search
[idxNN, dNN] =  rangesearch(pos, posquery, distlim);
maxNN = max(cellfun(@length, idxNN));
idxNN = cellfun(@(v) [v, zeros(1, maxNN-length(v))], idxNN, 'uni', 0);
idxNN = reshape([idxNN{:}], maxNN, [])';
dNN = cellfun(@(v) [v, nan(1, maxNN-length(v))], dNN, 'uni', 0);
dNN = reshape([dNN{:}], maxNN, [])';

% step distance
positives = [];
d = 5;
while d < distlim
    tempidx = idxNN.*(double(dNN<=d));
    nameNN = tempidx;
    nameNN(nameNN~=0) = idxName(nameNN(nameNN~=0));
    nNN = sum(logical(nameNN),2);
    nNNtemp = find(nNN>=length(uniName));
    nameNN = nameNN(nNNtemp,:);
    positive = false(size(nameNN,1),1);
    for i = 1:size(nameNN,1)
        uniNN = unique(nameNN(i,:));
        positive(i) = length(uniNN(uniNN~=0)) == length(uniName);
    end
    positive = nNNtemp(positive);
    
    for i = 1:length(positive)
        groups = tempidx(positive(i),:);
        groups(groups==0) = [];
        centroid = mean(pos(groups,:),1);
        positives = [positives; centroid, d, posquery(positive(i),:)];
    end
     d = d + 5;
end
[~, unipos] = unique(positives(:,4:5), 'rows');
positives = positives(unipos,:);

% visualization
if ~isempty(positives)
    hold on;
    D = unique(positives(:,3));
    for i = 1:length(D)
        tempidx = positives(positives(:,3)==D(i),:);
        for j = 1:size(tempidx,1)
            plot([tempidx(j,1)-D(i), tempidx(j,1)-D(i), tempidx(j,1)+D(i), tempidx(j,1)+D(i), tempidx(j,1)-D(i)],...
                [tempidx(j,2)+D(i), tempidx(j,2)-D(i), tempidx(j,2)-D(i), tempidx(j,2)+D(i), tempidx(j,2)+D(i)], 'y');
        end
    end

    for i = 1:size(positives,1)
        plot([positives(i,4)-distlim*2, positives(i,4)-distlim*2, positives(i,4)+distlim*2, positives(i,4)+distlim*2, positives(i,4)-distlim*2],...
            [positives(i,5)+distlim*2, positives(i,5)-distlim*2, positives(i,5)-distlim*2, positives(i,5)+distlim*2, positives(i,5)+distlim*2],...
            'w');
    end
else
    disp('No co-occurence of all reads within given distance.');
end

end