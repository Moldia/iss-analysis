%% 
SpotNormPrctile = 98;
nRounds = 5;

%% get intensity
spotColors = importdata('blobs.csv');
tilepos = spotColors.data(mod(spotColors.data(:,1),5)==1,[3, 9:10]);
spotColors = cat(3,...
    spotColors.data(mod(spotColors.data(:,1),5)==1,6:8),...
    spotColors.data(mod(spotColors.data(:,1),5)==2,6:8),...
    spotColors.data(mod(spotColors.data(:,1),5)==3,6:8),...
    spotColors.data(mod(spotColors.data(:,1),5)==4,6:8),...
    spotColors.data(mod(spotColors.data(:,1),5)==0,6:8));

%% global position
startpos = getcsvtilepos('Tiled.csv');
pos = zeros(length(tilepos),2);
for i = 1:length(startpos)
    pos(tilepos(:,1)==i,:) = tilepos(tilepos(:,1)==i,2:3) + startpos(i,2:3);
end

%% what used to correct crosstalk?
candidateGenes = {'KIAA0652' 'PTEN1' 'OXSM' 'TMEM8A' 'CCDC105'};
taglist = importdata('taglist_SangerMut.csv');
taglist = cellfun(@(v) strsplit(v, ','), taglist, 'UniformOutput', 0);
taglist = cat(1, taglist{:});
correctBarcodes = taglist(contains(taglist(:,2), candidateGenes),1);
candidateGenes = taglist(contains(taglist(:,2), candidateGenes),2);
codeNum = barcode2num(taglist(:,1)); 
[candidateGenes, correctBarcodes]

[~, candidateCodes] = barcode2num(correctBarcodes);
[~, rawCode] = max(spotColors, [], 2);
rawCode = squeeze(rawCode);
rawDecode = rawCode(:,1)*1e4 + rawCode(:,2)*1e3 + rawCode(:,3)*1e2 + rawCode(:,4)*1e1 + rawCode(:,5);
rawDecode = cellfun(@(v) find(v==codeNum), num2cell(rawDecode), 'UniformOutput', 0);
nnnn = cellfun(@isempty, rawDecode);
rawDecode(nnnn) = {0};
rawDecode = cell2mat(rawDecode);

forCorrection = ismember(rawCode, candidateCodes, 'rows');
nnz(forCorrection)

%% intensity clustering
spotColorsNorm = bsxfun(@rdivide, spotColors, prctile(spotColors, SpotNormPrctile));

BleedMatrix = zeros(3, 3, nRounds); % (Measured, Real, Round)
for r = 1:nRounds
    m = squeeze(spotColorsNorm(forCorrection,:,r)); % data: nCodes by nBases
    
    [Cluster, v, s2] = ScaledKMeans(m, eye(3));
    for i = 1:3
        BleedMatrix(:,i,r) = v(i,:) * sqrt(s2(i));
    end
end

figure;
for i = 1:nRounds
    subplot(2,3,i);
    imagesc(BleedMatrix(:,:,i));
    caxis([0 1]);
    title(sprintf('Round %d', i));
    set(gca, 'xtick', 1:4);
    set(gca, 'XTickLabel', {'A' 'C' 'G'});
    set(gca, 'ytick', 1:4);
    set(gca, 'yTickLabel', {'A' 'C' 'G'});
    if i==4
        xlabel('Actual')
        ylabel('Measured');
    end
end
subplot(2,3,6);
caxis([0 1]); 
axis off
colormap hot
colorbar

%% create code matrix
taglist = importdata('taglist_SangerMut.csv');
taglist = cellfun(@(v) strsplit(v, ','), taglist, 'UniformOutput', 0);
taglist = cat(1, taglist{:});

[~, code] = barcode2num(taglist(:,1));

idx = bsxfun(@plus, code, 0:3:3*nRounds-1);
idx = [reshape(idx', [], 1),...
    reshape(repmat(1:length(taglist),nRounds,1), [], 1)];

codemat = zeros(length(taglist), 3*nRounds);
idx = sub2ind(size(codemat), idx(:,2), idx(:,1));
codemat(idx) = 1;

%% code matrix with crosstalk
codematBled = zeros(size(codemat, 1), nRounds);

for i = 1:size(codemat,1)
    for r = 1:nRounds
        codematBled(i,3*(r-1)+(1:3)) = BleedMatrix(:,code(i,r),r);
    end
end

codematBledNorm = codematBled ./ sqrt(sum(codematBled.^2,2));

%% barcode matching
spotColorsFlat = spotColorsNorm(:,:);
spotIntensity = sqrt(sum(spotColorsFlat.^2,2));
spotColorsFlatNorm = bsxfun(@rdivide, spotColorsFlat, spotIntensity);
spotScores = spotColorsFlatNorm * codematBledNorm';

[spotScore, bestCode] = max(spotScores, [], 2);

figure,imagesc(spotScores);
set(gca, 'XTick', 1:length(taglist), 'XTickLabel', taglist(:,2), 'XTickLabelRotation', 90);

nnz(bestCode == rawDecode)

%% visualization
figure; set(gcf, 'Color', 'k', 'InvertHardcopy', 'off');

syms = symlist;

subplot(131);
plotall(taglist(rawDecode(rawDecode~=0),2), pos(rawDecode~=0,:), '');
update_legend(gca, taglist(:,2), syms(1:length(taglist)));
title('old decoding');

subplot(133);
plotall(taglist(bestCode(spotScore>0.8),2), pos(spotScore>0.8,:), '');
update_legend(gca, taglist(:,2), syms(1:length(taglist)));
title('new decoding');

subplot(132);
[name2, pos2] = getinsitudata('QT_0.4_0.001_details.csv');
plotall(name2(~strcmp(name2, 'NNNN')), pos2(~strcmp(name2, 'NNNN'),:), '');
update_legend(gca, taglist(:,2), syms(1:length(taglist)));
title('old decoding QT');

mkdir('individual');
for i = 1:length(taglist)
    subplot(131);
    update_legend(gca, taglist(i,2));
    legend off
    title(['old noQT ' taglist{i,2} ' n=' num2str(nnz(rawDecode(rawDecode~=0)==i))], 'Color', 'w');

    subplot(132);
    update_legend(gca, taglist(i,2));
    legend off
    title(['old QT=0.4 ' taglist{i,2} ' n=' num2str(nnz(strcmp(name2, taglist{i,2})))], 'Color', 'w');

    subplot(133);
    update_legend(gca, taglist(i,2));
    legend off
    title(['new ' taglist{i,2} ' n=' num2str(nnz(bestCode(spotScore>0.8)==i))], 'Color', 'w');

    saveas(gcf, ['individual\' taglist{i,2} '.png']);
end
