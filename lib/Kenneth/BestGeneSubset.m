function gl = BestGeneSubset(m0, gl);
% gl = BestGeneSubset(m0, gl);
% finds a set of genes that together identify cells to the best level of
% hierarhcy. m0 should be a MixNB structure, gl a list of genes to start
% with

m=m0;
Divisor = 50;
nIter = 100;
MinSumExp = 100; % a gene must have this many reads from all cells to be worth trying

[BinaryNames, ClassRoots] = BinarizeClassNames(m.ClassName);
PadNames = char(BinaryNames);
SimMat = zeros(m.nK);
for i=1:m.nK
    SimMat(i,i) = length(BinaryNames{i});
    for j=i+1:m.nK
        SimMat(i,j) = find(PadNames(i,:)~=PadNames(j,:), 1, 'first')-1;
        SimMat(j,i) = SimMat(i,j);
    end
end


DropThresh = 0;

AllMu = zeros(m.nG, m.nK);
for k=1:m.nK
    MyCells = (m.Class==k);
    AllMu(:,k) = (sum(m.x(:,MyCells),2) + m.RegN)/sum(MyCells + m.RegD);
end
% sx = m.x*m.w; % summed exp of all genes, nG by nK
% n = sum(m.w,1); % total cells in each class, 1 by nK 
% AllMu = bsxfun(@rdivide, sx+m.RegN, n+m.RegD)/Divisor; %nG by nK

if nargin<2
    CurrentSet = [];
else
    [~, CurrentSet] = ismember(gl, m.GeneName);
    CurrentSet = CurrentSet(:)'; % god i hate matlab
end

for i=1:nIter
    NewBest = CurrentSet;
    m = m0; % subsampled
    if ~isempty(Divisor) % subsample
        m.x = poissrnd(m0.x/Divisor);
    end

    OldScore = GeneSubsetPower(m, CurrentSet, SimMat, AllMu);

    % try adding
    WorthTrying = find(sum(m.x,2)>MinSumExp);
    Added = '';
    BestScore=OldScore;
    TryThese = setdiff(WorthTrying(:)', CurrentSet);
    Score = zeros(length(TryThese),1);
    parfor i=1:length(TryThese)
        gene = TryThese(i);
    %    for gene=setdiff(WorthTrying(:)', CurrentSet);
        TryGenes = [CurrentSet, gene];
        Score(i) = GeneSubsetPower(m, TryGenes, SimMat, AllMu);
%         if Score>BestScore
%             Added = gene;
%             Gain = Score-CurrentScore;
%             BestScore = Score;
%         end
%         if mod(gene,1000)==0
%             fprintf('.');
%         end
    end
    [BestScore, Best] = max(Score);
    if BestScore>OldScore
        Added = TryThese(Best);
        Gain = Score(Best) - OldScore;
    else
        Added = [];
    end
    
    if ~isempty(Added)
        CurrentSet = [CurrentSet, Added];
        fprintf(2, 'Added %s for gain of %f\n', m.GeneName{Added}, Gain);
            NewScore = GeneSubsetPower(m, CurrentSet, SimMat, AllMu);

    end
    
    % now try dropping
    
    GeneWorth = zeros(length(CurrentSet),1);
    parfor i=1:length(CurrentSet)
        gene=CurrentSet(i);
        TryGenes = setdiff(CurrentSet, gene);
        GeneWorth(i) = NewScore-GeneSubsetPower(m, TryGenes, SimMat, AllMu);
%         fprintf('Without %s: %f\n', m.GeneName{gene}, Score);
    end
    
    
% now print out some results, in a way athat can be cut-pasted into a
% spreadsheet

    fprintf('Gene\tScore\tMedianBestClass\n');
    for i=1:length(CurrentSet)
        g = CurrentSet(i);
        mu = zeros(m.nC,1);
        for c=1:m.nC
            mu(c) = median(m0.x(g,m.Class==c));
        end
        fprintf('%s \t%.4f \t%.1f\n', m.GeneName{g}, GeneWorth(i), max(mu));
    end

%     for i=1:length(CurrentSet)
%         fprintf('%s\t%f\n', m.GeneName{CurrentSet(i)}, GeneWorth(i));
%     end
    % now sort into order
    [WorthsSorted, order] = sort(GeneWorth, 'descend');
    
    if min(GeneWorth)<DropThresh;
        fprintf(2, '\nDropping %s for gain of %f\n', m.GeneName{CurrentSet(order(end))},...
            WorthsSorted(end));
        CurrentSet(order(end))=[];
    end

        
    gl = m.GeneName(CurrentSet);

    fprintf('\nTotal Score %f\n\n', NewScore)

end


return
%%
for g = gl(:)'
    Max = max(gBoth.Exp(g));
    figure(101); gBoth.SortByClass.BoxPlot(g); ylim([0 Max]);
    figure(102); gPCs0.SortByClass.BoxPlot(g); ylim([0 Max]);
    figure(103); gINs0.SortByClass.BoxPlot(g); ylim([0 Max]);
    pause
end

