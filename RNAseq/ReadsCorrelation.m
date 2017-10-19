%% different lanes

FPKM = cell(4,4);

for s = 25:28
    for l = 1:4
        data = importdata(fullfile('UngappedMapping\fpkm_cufflinks',...
            ['S' num2str(s) '_L' paddigits(l, 3)],...
            'fpkm_sipmle.csv'));
        
        name = data.textdata(:,2);
        fpkm = data.data;
        
        FPKM{s-24,l} = [name, num2cell(fpkm)];
                
    end
end

save UngappedMapping\fpkm_cufflinks\FPKM.mat FPKM

%% find the same genes
allGenes = cellfun(@(v) v(:,1), FPKM, 'uni',0);
allGenes = unique(cat(1, allGenes{:}));
allFPKM = cell(4,4);
meanFPKM = [];
stdFPKM = [];

for s = 1:4
    lanes = zeros(numel(allGenes), 4);
    for l = 1:4
        order = cellfun(@(v) find(strcmp(v, allGenes)), FPKM{s,l}(:,1));
        lanes(order,l) = cell2mat(FPKM{s,l}(:,2));
    end
    allFPKM{s,l} = lanes;
    meanFPKM(:,s) = mean(lanes,2);
    stdFPKM(:,s) = std(lanes, 0, 2);
end

figure, hold on;
bar(meanFPKM(:,1));
e = errorbar(meanFPKM(:,1), stdFPKM(:,1), '.');
e.CapSize = 5;
        


