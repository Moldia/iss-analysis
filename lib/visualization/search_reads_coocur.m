function search_reads_coocur(name, pos, distlim, query)
% find where reads occur together with given query gene
% Xiaoyan, 2017


[uniName, ~, idxName] = unique(name);
if nargin > 3
    idxQ = find(strcmp(uniName, query));
else
    countName = hist(idxName, 1:length(uniName));
    [~, idxQ] = sort(countName, 'ascend');
    idxQ = idxQ(1);
end


positives = [];
d = 5;
nNNtypes = 0;
while d < distlim
    posquery = pos(idxName==idxQ,:);
    [idxNN, ~] =  rangesearch(pos, posquery, d);
    nameNN = cellfun(@(v) unique(idxName(v)), idxNN, 'uni', 0);
    nNNtypes = cellfun(@length, nameNN);
    wNN = nNNtypes > 1;
    centroid = cellfun(@(v) mean(pos(v,:)), idxNN(wNN), 'uni', 0);
    positives = [positives;...
        (reshape([centroid{:}], 2, []))', repmat(d, length(centroid), 1),...
        nNNtypes(wNN), posquery(wNN,:)];
    d = d + 5;
end
[~, unipos] = unique(positives(:,4:6), 'rows');
positives = positives(unipos,:);

% visualization
hold on;
D = unique(positives(:,3));
for i = 1:length(D)
    temp = positives(positives(:,3)==D(i),:);
    for j = 1:size(temp,1)
        plot([temp(j,1)-D(i), temp(j,1)-D(i), temp(j,1)+D(i), temp(j,1)+D(i), temp(j,1)-D(i)],...
            [temp(j,2)+D(i), temp(j,2)-D(i), temp(j,2)-D(i), temp(j,2)+D(i), temp(j,2)+D(i)], 'y');
    end
end

for i = 1:size(positives,1)
    plot([positives(i,5)-distlim*2, positives(i,5)-distlim*2, positives(i,5)+distlim*2, positives(i,5)+distlim*2, positives(i,5)-distlim*2],...
        [positives(i,6)+distlim*2, positives(i,6)-distlim*2, positives(i,6)-distlim*2, positives(i,6)+distlim*2, positives(i,6)+distlim*2],...
        'w');
end