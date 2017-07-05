function [Genes, Codes, MaxScore, Intensity] = iss_call_spots(o, SpotColors0, IsolatedSpots)
% calls spots to codes for in-situ sequencing. Run this after generating
% SpotColors and IsolatedSpots with InSituFindSpots.m. CodeFile is a
% filename pointing to the "codebook_unique.csv" file.
% 
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
  
nChans = o.nBP+2;

SpotColors = bsxfun(@rdivide, SpotColors0, prctile(SpotColors0, o.SpotNorm(1), o.SpotNorm(2)));


% now we cluster the intensity vectors to estimate the Bleed Matrix
BleedMatrix = zeros(o.nBP,o.nBP,o.nRounds); % (Measured, Real, Round)
for r =1:o.nRounds
    m = sq(SpotColors(IsolatedSpots,:,r)); % data: nCodes by nBases
    
    [Cluster, v, s2] = ScaledKMeans(m, eye(4));
    for i=1:4
        BleedMatrix(:,i,r) = v(i,:) * sqrt(s2(i));
    end
end

if o.Graphics
    figure(98043765); clf
    for i=1:5
        subplot(2,3,i); 
        imagesc(BleedMatrix(:,:,i)); 
        caxis([0 1]); 
        title(sprintf('Cycle %d', i)); 
        set(gca, 'xtick', 1:4);
        set(gca, 'XTickLabel', o.bpLabels);
        set(gca, 'ytick', 1:4);
        set(gca, 'yTickLabel', o.bpLabels);
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
    save BleedMatrix BleedMatrix
end

% now load in the code book and apply bleeds to it
codebook_raw = importdata(o.CodeFile);
CharCode = codebook_raw.textdata(2:end,5);
GeneName = codebook_raw.textdata(2:end,3);
nCodes = size(CharCode,1)-2; % bit of a hack to get rid of Sst and Npy

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r = 1:o.nRounds
    NumericalCode(:,r) = codebook_raw.data(1:nCodes,(r-1)*nChans + (3:nChans))*(1:o.nBP)';
end

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:o.nRounds
        BledCodes(i,(1:o.nBP) + (r-1)*o.nBP) = BleedMatrix(:, NumericalCode(i,r), r);
        UnbledCodes(i,NumericalCode(i,r) + (r-1)*o.nBP) = 1;
    end
end

NormBledCodes = bsxfun(@rdivide, BledCodes, sqrt(sum(BledCodes.^2,2)));
FlatSpotColors = SpotColors(:,:);
Intensity = sqrt(sum(FlatSpotColors.^2,2));
NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, Intensity);
SpotScores = NormFlatSpotColors * NormBledCodes';

[MaxScore, BestCode] = max(SpotScores,[],2);
Codes = CharCode(BestCode);
Genes = GeneName(BestCode);

