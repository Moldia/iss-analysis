function [pCellClass, ClassNames, pSpotCellSparse] = iss_call_cells(o, SpotYX, SpotGenes, gSet, CellMap)
% [pCellClass, pSpotCell] = iss_call_cells(o, SpotYX, SpotGenes, gSet, CellMap)
%  
% Cell calling via negative binomial model
% 
% inputs: o, options structure
% SpotYX: coordinates of RNA detections (nSpots*2, y then x coord)
% SpotGenes: genes of each detected spot (nSpots*1 cell array of strings)
% gSet: GeneSet structure containing results of scRNAseq clustering
% CellMap: output of Dapi segmentation (
%
% outputs: 
% pCellClass: posterior probability of each cell to be in each class (nCells*nClasses)
% pSpotCell: posterior probability of each spot to be in top 5 neighboring cells (nSpots * nCells, sparse)
% note that the last class is a zero-expressing cell; the last cell is background


%
% if ~exist('gSet')
%     load gSet
%     load ../Jens/Slc32a1CA1/gAll4
%     load ../Jens/Slc32a1CA1/gFinal.mat
%     ddd = importdata('../Jens/Slc32a1CA1/Slc32a1 classes all.xlsx');
%     OldINNames = deblank(ddd.textdata(2:end,1));
%     NewINNames = deblank(ddd.textdata(2:end,2));
%     gSet0 = ReClass([gPCs gNonNeuron g], containers.Map(OldINNames, NewINNames));
%     gSet = gSet0.CellSubset('', 'Contaminated');
% end
% gSet=hAll.CellSubset('', 'Unclassified');

% SpotYX = FinalYX(FinalMaxScore>.9,:);
% SpotGenes = FinalGenes(FinalMaxScore>.9);
% SpotGenes(strcmp('Lphn2',SpotGenes))={'Adgrl2'};
%o.CellSize = 15; % cell radius in Pixels, for Gaussian penalty
rp = regionprops(CellMap);
CellYX = fliplr(vertcat(rp.Centroid)); % convert XY to YX
CellArea0 = vertcat(rp.Area); 
MeanCellArea = mean(CellArea0);
%RelCellArea = [CellArea0/MeanCellArea; 1]; 
RelCellArea = [ones(size(CellArea0)) ; 1]; % last "cell" is background signal 

o.nNeighbors = 3; % how many neighboring cells to consider for each spot
o.rSpot = 2; % gamma dist prior for each spot
o.rGene = 2; % gamma dist shape for scaling whole gene
o.MisreadDensity = 1e-4; % dodgy reads per pixel
o.Inefficiency = .2; % how much to scale scRNAseq down by to match iss
o.SpotReg = .2; % how much to add to expression so as not to get infities
o.InsideCellBonus = 3; % additional log likelihood for a gene inside the segmented region

o.CellCallTolerance = .01; % converges when no probabilities have changed more than this


o.CellCallShowCenter = [1670 250];
o.CellCallShowRad = 200;
o.ExampleCellCenter = [1670 250];

%o.ClassPrior = [repmat(.7/49,1,49) .1 .1 .1]; % last is PC, nonN, Zeros

%% get arrays ready

% SpotGene(nS): which gene is each spot
% MeanClassExp(nK,nG): mean expression of each gene in each class
% Neighbors(nS, nN): closest neighboring cells for each spot
% D(nS, nN): distance penalty for each of these
% GeneNames(nG): name of each gene
% ClassNames(nK): name of each class

[GeneNames, ~, SpotGene] = unique(SpotGenes);
TotGeneSpots = accumarray(SpotGene,1);
[ClassNames] = vertcat(unique(gSet.Class, 'stable'), {'Zero'});

nG = length(GeneNames);
nK = length(ClassNames); % last is zero-expression
nC = size(CellYX,1)+1; % last is misreads
nS = size(SpotYX,1);
nN = o.nNeighbors+1; % last is misreads (always a neighbor)

o.ClassPrior = [.5*ones(1,nK-1)/nK .5];

ClassDisplayNames = ClassNames;
% ClassDisplayNames(1:5) = {'PC'};
% ClassDisplayNames(6:22) = {'NonNeuron'};
%%
MeanClassExp = zeros(nK, nG);
gSub = gSet.GeneSubset(GeneNames);
for k=1:nK-1 % don't include last since it is zero-expression class
    MeanClassExp(k,:) = o.Inefficiency * mean(gSub.ScaleCell(1).CellSubset(ClassNames{k}).GeneExp,2)';
end
MeanClassExp = MeanClassExp + o.SpotReg;
lMeanClassExp = log(MeanClassExp); 

% now find each spot's neighboring cells and distances
[Neighbors, Dist] = knnsearch(CellYX, SpotYX, 'K', nN);
Neighbors(:,end) = nC; % set last neighbor to misreads
%D = -Dist.^2/(2*o.CellSize.^2) - log(2*pi*o.CellSize.^2);
D = -Dist.^2./(2*RelCellArea(Neighbors)*MeanCellArea/pi) - log(2*MeanCellArea); % don't normalize: bigger cells express more
D(:,end) = log(o.MisreadDensity); % this is log likelihood of misread

% any inside cell radius given a bonus
SpotInCell = IndexArrayNan(CellMap, SpotYX');
if Neighbors(SpotInCell>0,1)~=SpotInCell(SpotInCell>0)
    error('a spot is in a cell not closest neighbor!');
end
D(SpotInCell>0, 1) = D(SpotInCell>0, 1) + o.InsideCellBonus;



% LogClassPrior = zeros(1,nK);
% LogClassPrior(1:nK-1) = log((1-o.ZeroPrior)/(nK-1));
LogClassPrior = log(o.ClassPrior);



% Dist0 = zeros(nS, nC);
% for c=1:nC
%     dt = bwdist(CellMap==c);
%     Dist0(:,nC) = dt(sub2ind(size(dt), SpotYX(:,1), SpotYX(:,2)));
% end
% [sDist, oDist] = sort(Dist0,2, 'ascend');
% Dist = zeros(nS, nN);
% Dist(:,1:nN-1) = sDist(:,1:nN-1);
% Neighbors = oDist(:,1:nN-1);
%% variables for main loop

pSpotCell = sparse(nS, nN); % prob each spot goes to each neighboring cell: last assigned to noise
pCellClass = zeros(nC, nK); % prob each cell goes to each class: last has zero expression
% initialize arrays

% each spot goes to closest neighbor
pSpotCell(:,1)=1;
% gammas start of as priors

% SpotGammaShape = zeros(nC, nK, nG); % gamma dist parameters for each spot
% SpotGammaRate = zeros(nC, nK, nG); % gamma dist parameters for each spot
SpotGammaShape(:,:,:) = o.rSpot;
SpotGammaRate(:,:,:) = o.rSpot;
% eSpotGamma = SpotGammaShape./SpotGammaRate; % epectation of gamma (nC,nK,nG)
% elSpotGamma = psi(SpotGammaShape) - log(SpotGammaRate); % expectation of log gamma

eSpotGamma = ones(nC, nK, nG);
elSpotGamma = ones(nC, nK, nG)*psi(1); % start with r=1 prior, no spots

eGeneGamma = ones(nG,1); % start with just 1

pSpotCellOld = zeros(nS, nN);
for i=1:100
    % CellGeneCount(nC, nG): number of copies of each gene in each cell
    CellGeneCount = zeros(nC,nG);
    for n=1:nN-1
        c = Neighbors(:,n);
        CellGeneCount = CellGeneCount + accumarray([c, SpotGene], pSpotCell(:,n), [nC,nG]);
    end

    %% call cells
    % wCellClass(nC, nK)
    % p(nC, nK, nG)
    ScaledExp = reshape(MeanClassExp,[1 nK nG]) .* reshape(eGeneGamma,[1 1 nG]) .* RelCellArea;
    p = ScaledExp ./ (o.rSpot + ScaledExp);
    %p = 1./(1+o.rSpot./bsxfun(@times, MeanClassExp, eGeneGamma')); % p(nK, nG)
    
    %wCellClass(nC,nK)
    wCellClass = sum(reshape(CellGeneCount,[nC 1 nG]).*log(p) + o.rSpot*log(1-p),3) + LogClassPrior;
    
    %pCellClass(nC, nK)
    pCellClass = LogLtoP(wCellClass')';

%     if o.Graphics
%         figure(3985470); cla
%         iss_make_figure(o, FinalYX(ShowMe,:), FinalGenes(ShowMe)); 
%         hold on
%         scatter(CellYX(:,2), CellYX(:,1), 100, 'w');
%         scatter(CellYX(:,2), CellYX(:,1), 25, 'w');
%         [~, BestClass] = max(pCellClass(1:end-1,:),[],2);
%         text(CellYX(:,2), CellYX(:,1), ClassNames(BestClass), 'color', 'w');
%     end
    
    %% call gammas
%     SpotGammaShape = o.rSpot + bsxfun(@times, reshape(CellGeneCount,[nC,1,nG]), reshape(pCellClass,[nC,nK,1]));
%     SpotGammaRate = o.rSpot + bsxfun(@times, reshape(MeanClassExp,[1,nK,nG]), reshape(pCellClass,[nC,nK,1]));
    
    %eSpotGamma(nC, nK, nG);
    ScaledMean = RelCellArea.*reshape(MeanClassExp,[1 nK nG]);
    eSpotGamma = (o.rSpot+reshape(CellGeneCount,[nC 1 nG]))./(o.rSpot + ScaledMean);
    elSpotGamma = psi(o.rSpot+reshape(CellGeneCount,[nC 1 nG])) - log(o.rSpot + ScaledMean); % expectation of log gamma

%    eGeneGamma = (o.rGene + TotGeneSpots)./(o.rGene + sum(pCellClass*MeanClassExp,1)');
    % to count non-background expression of each gene
    BackgroundSpots = accumarray(SpotGene, pSpotCell(:,end), [nG 1]);
    % total predicted by other models
    TotPredicted = sum(shiftdim(sum(eSpotGamma.*pCellClass.*RelCellArea,1),1).*MeanClassExp,1)';
    eGeneGamma = (o.rGene + TotGeneSpots - BackgroundSpots)./(o.rGene + TotPredicted);
    if 0
        for gg=1:nG; fprintf('%s:\t%f\n', GeneNames{gg}, eGeneGamma(gg)); end
    end
    
    %% call spots
    % wSpotCell(nS, nN)
    aSpotCell = zeros(nS, nN);
    for n=1:nN-1 % don't include misread possibility
        c = Neighbors(:,n);
        aSpotCell(:,n) = sum(pCellClass(c,:) .* lMeanClassExp(:,SpotGene)',2) + ...
            sum(pCellClass(c,:) .* bi(elSpotGamma, c, 1:nK, SpotGene), 2);
    end
    wSpotCell = aSpotCell + D ;
    
    pSpotCell = LogLtoP(wSpotCell')';
    MeanProbChanged = max(abs(pSpotCell(:)-pSpotCellOld(:)));
    fprintf('Iteration %d, mean prob change %f\n', i, MeanProbChanged)
    if MeanProbChanged<o.CellCallTolerance
        break
    end
    pSpotCellOld = pSpotCell;
    
end

%% make dense array output

pSpotCellSparse = sparse(repmat(1:nS,1,nN)', Neighbors(:), pSpotCell(:));

%% diagnostics
if ~isempty(o.CellCallShowCenter)
    figure(3985471); cla
    iss_make_figure(o, SpotYX, SpotGenes); 
    hold on
    [~, BestNeighb] = max(pSpotCell,[],2);
    %BestNeighb = ones(nS,1);
    SpotBestNeighb = bi(Neighbors,(1:nS)',BestNeighb(:));
    rn = SpotBestNeighb<nC & sum(abs(SpotYX-o.CellCallShowCenter).^2,2)<o.CellCallShowRad.^2;
    plot([SpotYX(rn,2) , CellYX(SpotBestNeighb(rn),2)]', ...
        [SpotYX(rn,1) , CellYX(SpotBestNeighb(rn),1)]', 'Color', [.3 .3 .3]);
    [~, BestClass] = max(pCellClass(1:end-1,:),[],2);            
    text(CellYX(:,2), CellYX(:,1), ClassDisplayNames(BestClass), 'color', 'r', 'fontsize', 6);
%     input('press key', 's');
    axis(reshape(fliplr([o.CellCallShowCenter;o.CellCallShowCenter])+[-o.CellCallShowRad; o.CellCallShowRad], 1, 4));

end


%%
% diagnostics:
if ~isempty(o.ExampleCellCenter)
    [~, MyCell] = min(sum((CellYX-o.ExampleCellCenter).^2,2));
    fprintf('------------------ Cell %d: -----------------\n', MyCell);
    for ss=find((Neighbors(:,1)==MyCell))'
        fprintf('Spot %d: %s, with prob %f\n', ss, GeneNames{SpotGene(ss)}, pSpotCell(ss,1));
    end
    fprintf('-- Total Gene Count --\n');
    for gg=find(CellGeneCount(MyCell,:)>1e-3) 
        fprintf('%s:\t%f\n', GeneNames{gg}, CellGeneCount(MyCell,gg)); 
    end
    fprintf('-- Class Posteriors --\n');
    for cc=find(pCellClass(MyCell,:)>1e-3)
        fprintf('%s:\t%e\n', ClassDisplayNames{cc}, pCellClass(MyCell,cc)); 
    end
end
%% 
%%
return
% figure(398547); clf
% iss_make_figure(o, FinalYX(ShowMe,:), FinalGenes(ShowMe)); hold on
% [~, BestClass] = max(pCellClass(1:end-1,:),[],2);
% text(CellYX(:,2), CellYX(:,1), ClassNames(BestClass), 'color', 'w');
% % Colors = [HsvNotYellow(nK-4); .5 .5 .5 ; .8 .8 .8 ; .1 .1 .1 ; 0 0 0];
% scatter(CellYX(:,2), CellYX(:,1), 75, Colors(BestClass(1:end-1),:));
% scatter(CellYX(:,2), CellYX(:,1), 50, Colors(BestClass(1:end-1),:));
% set(gca, 'Color', 'k');