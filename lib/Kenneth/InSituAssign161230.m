% T = readtable('../../InSituSequencing/HC Cortex cortex probes only/CellBlobs_160408_1610222UCL_hippo_1-GCaMP.csv');
% 
% Fig = openfig('../../InSituSequencing/HC Cortex cortex probes only/Plotted.fig');
% ChangeGeneSymbols;
% ShowGenes('', 'Slc1a2');
% Fig2 = openfig('../../InSituSequencing/HC Cortex cortex probes only/Plotted.fig');
% ChangeGeneSymbols;
% ShowGenes('', 'Slc1a2');
% 
% xr = [980 4038];
% yr = [7644 10100];
% load ../Jens/Htr3aCA1/gINs161023
% load ../Jens/Htr3aCA1/mNonINs
%%
T = readtable('PerCell/PerCellCount_ROI1_NoCell0NoEmptyCells.csv');

ZoomFig = openfig('161230WithSstNpy.fig');
ChangeGeneSymbols;
ShowGenes('');
%camroll(90);

% BoxFig = openfig('161230WithSstNpy.fig');
% %figure;%openfig('../../InSituSequencing/GCamp HC Cortex all probes/Plotted_CellOutlines.fig');
% ChangeGeneSymbols;
% ShowGenes('');
% camroll(90);
% legend('off');

PieFig = figure;
set(PieFig, 'color', [0 0 0]);
%camroll(90);
axis off

% top left
xr = [4082   6063];
yr = [1564   4153];

% bottom right
xr = [1847   4729];
yr = [9400  12019];

% big bottom right
xr = [2382 6461];
yr = [7176 12060];

% everywhere
xr = [2020 6200];
yr = [1214 12000];

BoxSize = 50;
load ../Jens/Htr3aCA1/gINs161023
load ../Jens/Htr3aCA1/mNonINs

%%
% colors
non_neuron = hsv2rgb([0 0 1]);
pc_or_in =   hsv2rgb([.4 .5 .5]);
pc =        hsv2rgb([1/3 1 1]);
pc2 =       hsv2rgb([1/4 1 .9]);
in_general = hsv2rgb([2/3 1 1]);

sst =   hsv2rgb([.55 1 1]);
pvalb = hsv2rgb([.7 .8 1]);
ngf =   hsv2rgb([.85 1 1]);
cnr1 =  hsv2rgb([ 1 1 1]);
vip =   hsv2rgb([ .15 1 1]);
cxcl14= hsv2rgb([.1 1 .6]);

proj = hsv2rgb([1/12 1 1]);
NickNamesColors = {...
'Pvalb Basket', pvalb; ...
'Sst Bistrat', sst; ...
'Sst SeptProj', sst; ...
'Sst GIN', sst; ...
'Sst OLM', sst; ...
'Bcl11b Ntng1 Ndnf', ngf; ...
'Chrm2 Trilaminar', sst; ...
'Pvalb Chandelier', pvalb; ...
'Ngf MGE/Ivy', ngf; ...
'Ngf CGE', ngf; ...
'Ngf CGE Npy++', ngf; ...
'Ngf Cxcl14',ngf;  ...
'Cck Vip Cxcl14', cnr1; ...
'IS2', vip; ...
'Cck Vip ', cnr1; ...
'Cck Vip Tac2', cnr1; ...
'Cck Calca', cnr1; ...
'Cck Npy Tac2', cxcl14; ...
'Cck Npy Slc17a8', cxcl14; ...
'Cck Npy Cpne5', cxcl14; ...
'Cck Calb1', cxcl14; ...
'Cck Lypd1', cxcl14; ...
'Cck Lxn', cxcl14; ...
'IS1 Sln', vip; ...
'IS1 Tac2', vip; ...
'IS1 Igfbp6', vip; ...
'IS3 Igfbp4', vip; ...
'IS3 Crh', vip; ...
'IS3 Myl1', vip; ...
'PC CA2', [0 .3 0]; ...
'PC CA1', pc; ...
'PC CA1', pc; ...
'Endo', non_neuron; ... % CA2 in a different color
'Astro', non_neuron; ...
'Oligo', non_neuron...
};

NickNames=NickNamesColors(:,1);
Colors = NickNamesColors(:,2);
rand('state', 1);
for i=1:length(Colors)
    Colors{i} = max(min(Colors{i}*.8 + .1 + (rand(1,3)-.5)*0.6,1),0);
end

d1 = .02;
d2 = .4;
dy = .1;
nPerCol = 4;
nCols = floor(length(Colors)/nPerCol);
figure(237890); clf; hold on
ReorderLegend = [1 8 2:5 7 6 9:13 15:26 14 27:35];
for j=0:length(Colors)-1
    i = ReorderLegend(j+1);
%     plot(d2*floor(j/18), -mod(j, 18)-floor(j/18), 'o', 'markeredgecolor', 'none', 'markerfacecolor', Colors{i});
%     text(d1 + d2*floor(j/18), -mod(j, 18)-floor(j/18), NickNames{i});
    plot(d2*floor(j/nPerCol), -dy*mod(j, nPerCol), 'o', 'markeredgecolor', 'none', 'markerfacecolor', Colors{i});
    h = text(d1 + d2*floor(j/nPerCol), -dy*mod(j, nPerCol), NickNames{i}); %-floor(j/nPerCol)
    set(h, 'fontsize', 15);

end
set(gcf, 'Color', [1 1 1]*.8);
set(gca, 'Color', [1 1 1]*.8);
xlim([0 d2*nCols+.8])
axis off
DoubleColFig(nPerCol,8*nCols);    
set(gcf, 'InvertHardcopy', 'off');
%print -dpng -r300 ../../Talks/sfn' poster 2016'/Cell_calling_legend.png



mAll = struct;
mAll.x = horzcat(mBest.x, mNonINs.x);
mAll.nC = mBest.nC + mNonINs.nC;
mAll.nK = mBest.nK + mNonINs.nK;
mAll.Class = vertcat(mBest.Class, mNonINs.Class + mBest.nK);
mAll.GeneName = mBest.GeneName;
mAll.OptimalClassOrder = vertcat(mBest.OptimalClassOrder, mNonINs.OptimalClassOrder);
mAll.ClassName = vertcat(mBest.ClassName, mNonINs.ClassName);
mAll.r = mBest.r;

%%

nG0 = length(T.Properties.VariableNames)-3;
GeneNames0 = T.Properties.VariableNames(2:1+nG0);

%KeepGenes = find(~ismember(GeneNames0, {'Slc1a2', 'Plp1', 'NNNN', 'Lphn2'}));
% KeepGenes = find(~ismember(GeneNames0, {'Slc1a2', 'NNNN', 'Lphn2'}));
GeneNames = GeneNames0;
GeneNames{1} = '3110035E14Rik';
GeneNames{39} = 'Adgrl2';
nG = length(GeneNames);

gx0 = table2array(T);
gx = gx0(:,2:end-2);

nSpots = sum(gx,2);
xPos = gx0(:,end-1);
yPos = gx0(:,end);

[~, order] = sort(nSpots, 'descend');
%% make mean expression of selected genes
    MyGenes = cell(nG,2);
    MyGenes(:,1) = GeneNames;
    MyGenes(:,2) = num2cell(gx(i,:));

GeneNames = MyGenes(:,1);
[~, GeneIDs] = ismember(GeneNames, mAll.GeneName);
nMy = length(GeneIDs);

MyMeans = zeros(nMy,mAll.nK);

RegN = 10; RegD = 10;
for i=1:mAll.nK
    k = strmatch(mAll.OptimalClassOrder(i), mAll.ClassName);
    MyMeans(:,i) = (RegN+sum(mAll.x(GeneIDs, mAll.Class==k),2))/(RegD + sum(mAll.Class==k));
end

%% now do it
BoxSize=10;

CellIDs = 31:length(order); %[24 105 176 199 543 1028 1097 1102 1105 1109];

figure(ZoomFig)
%DoubleColFig(4,4);
ch = get(gca, 'children');
%set(ch(1:end-1), 'markersize', 8);

Interactive = 1;
% figure(BoxFig); xlim(xr); ylim(yr);
% figure(PieFig); xlim(xr); ylim(yr);
for ii=CellIDs; %1:length(order)
    i = order(ii);
    if nSpots(i)<=2; break; end
%     if ~(xPos(i)>xr(1) & xPos(i)<xr(2) & yPos(i)>yr(1) & yPos(i)<yr(2))
%         continue;
%     end
    fprintf('Cell at %.0f, %.0f: ', xPos(i), yPos(i));
    MyGenes = cell(nG,2);
    MyGenes(:,1) = GeneNames;
    MyGenes(:,2) = num2cell(gx(i,:));
    for j=1:nG
        if MyGenes{j,2}>0
            fprintf('%s %d ', GeneNames{j}, MyGenes{j,2});
        end
    end
    fprintf('\n');
    

    drawnow
    if Interactive
        x0 = (xPos(i)+22707)/4;
        y0 = (yPos(i)+6097)/4; 
        Scale = .25;
        figure(ZoomFig); xlim([-1 1]*BoxSize+x0); ylim([-1 1]*BoxSize+y0)
        %figure(BoxFig); BoxHandle(i) = plot([-1 1 1 -1 -1]*BoxSize+xPos(i), [-1 -1 1 1 -1]*BoxSize+yPos(i), 'w', 'linewidth', 2);
        figure(PieFig); CellCalla(MyMeans, gx(i,:), NickNames, Colors); set(gcf, 'color', 'w')

        inp = input('p to print to file', 's');
        if strcmp(inp,'p')
            figure(ZoomFig);
            print(sprintf('pix/Zoom %d', ii), '-dpng', '-r300');
            figure(BoxFig);
            print( sprintf('pix/Box %d', ii), '-dpng', '-r300');
            figure(PieFig);
            print(sprintf('pix/Pie %d', ii), '-dpng', '-r300');
        end
%         figure(BoxFig); delete(BoxHandle(i));
    else
        figure(PieFig);
        PieHandle{ii} = CellCalla(MyMeans, gx(i,:), NickNames, Colors, 20, [xPos(i), yPos(i)]);

    end

end