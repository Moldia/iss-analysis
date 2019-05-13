data = importdata('BlastOverview.csv');
q = data.textdata(2:end,1);
s = data.textdata(2:end,2);

figure;

%% self-dimers
dimer_self = false(length(q),1);
for i = 1:length(q)
    dimer_self(i) = strcmp(q{i},s{i});
end
dimer_hetero = find(~dimer_self);
dimer_self = find(dimer_self);

DimerSelf = zeros(max(reshape(data.data(:,5:8),[],1)));
for i = 1:length(dimer_self)
    DimerSelf(data.data(dimer_self(i),5):data.data(dimer_self(i),6),...
        data.data(dimer_self(i),8):data.data(dimer_self(i),7)) = ...
        DimerSelf(data.data(dimer_self(i),5):data.data(dimer_self(i),6),...
        data.data(dimer_self(i),8):data.data(dimer_self(i),7)) + 1;
end

subplot(221);
bh = bar3(DimerSelf, 1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
    bh(i).EdgeAlpha = 0.2;
end
view(2)
axis normal
colormap hot
colorbar
title('self dimers, base position')

%% heterodimers
DimerHetero = zeros(max(reshape(data.data(:,5:8),[],1)));
for i = 1:length(dimer_hetero)
    DimerHetero(data.data(dimer_hetero(i),5):data.data(dimer_hetero(i),6),...
        data.data(dimer_hetero(i),8):data.data(dimer_hetero(i),7)) = ...
        DimerHetero(data.data(dimer_hetero(i),5):data.data(dimer_hetero(i),6),...
        data.data(dimer_hetero(i),8):data.data(dimer_hetero(i),7)) + 1;
end

subplot(222);
bh = bar3(DimerHetero, 1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
    bh(i).EdgeAlpha = 0.2;
end
view(2)
axis normal
colormap hot
colorbar
title('hetero-dimers, base position')

%% heterodimer probe pairs
Pairs = [q(dimer_hetero),s(dimer_hetero)];
[name_uni,~,idx_name] = unique(Pairs(:));
name_uni = cellfun(@(v) strrep(v,'_','\_'),name_uni,'uni', 0);
DimerPairs = hist3(reshape(idx_name,[],2), [{1:length(name_uni)},{1:length(name_uni)}]); 

subplot(223);

L = linkage(DimerPairs,'average','cityblock');
dendrogram(L,0,'Orientation','right')
idx_probe = get(gca,'YTickLabel');
idx_probe = str2num(idx_probe);
idx_probe = flipud(idx_probe);

bh = bar3(DimerPairs(idx_probe,idx_probe), 1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
    bh(i).EdgeAlpha = 0.2;
end
view(2)
axis normal
colormap hot
colorbar
title('hetero-dimers, probe names')
set(gca,'XTick',1:length(name_uni),'XTickLabel',name_uni(idx_probe),...
    'YTick',1:length(name_uni),'YTickLabel',name_uni(idx_probe),...
    'XTickLabelRotation',90);

%% heterodimer, position-probe
DimerHeteroPos = zeros(length(dimer_hetero),max(reshape(data.data(:,5:8),[],1)));
for i = 1:length(dimer_hetero)
    DimerHeteroPos(i,data.data(dimer_hetero(i),5):data.data(dimer_hetero(i),6)) = ...
        DimerHeteroPos(i,data.data(dimer_hetero(i),5):data.data(dimer_hetero(i),6)) + 1;
end
[name_uni,~,idx_name] = unique(q(dimer_hetero));
name_uni = cellfun(@(v) strrep(v,'_','\_'),name_uni,'uni', 0);

DimerSpec = zeros(length(name_uni),max(reshape(data.data(:,5:8),[],1)));
for i = 1:length(name_uni)
    DimerSpec(i,:) = sum(DimerHeteroPos(idx_name==i,:),1);
end
DimerSpec = log2(DimerSpec+1);

subplot(224);

L = linkage(DimerSpec','weighted');
dendrogram(L,0,'Orientation','right')
idx_pos = get(gca,'YTickLabel');
idx_pos = str2num(idx_pos);
idx_pos = flipud(idx_pos);

L = linkage(DimerSpec(:,idx_pos),'weighted');
dendrogram(L,0,'Orientation','right')
idx_probe = get(gca,'YTickLabel');
idx_probe = str2num(idx_probe);
idx_probe = flipud(idx_probe);

bh = bar3(DimerSpec(idx_probe,idx_pos), 1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
    bh(i).EdgeAlpha = 0;
end
view(2)
axis normal
colormap hot
hc = colorbar;
set(hc,'Ytick',0:5,'YTicklabel',2.^(0:5));
title('hetero-dimers, position - probe names')
set(gca,'YTick',1:size(DimerSpec,1),'YTickLabel',name_uni(idx_probe),...
    'XTick',1:size(DimerSpec,2),'XTickLabel',idx_pos,...
    'XTickLabelRotation',90);
