%% top level script for analyzing in situ sequencing data
% specialized for dataset 170315 at 20x
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
% set parameters to default values
o = iss_options;
o.InputDirectory = '';
o.OutputDirectory = '';
o.TileDirectory = 'F:\In process\170315_161220KI_4-3\Preprocessing\';
o.CodeFile = '..\codebook_unique_full.csv';
o.bpLabels = {'A', 'C', 'G', 'T'};
o.TileFiles = cell(5,1,1);
o.TileFiles{1,1,1} = 'F:\In process\170315_161220KI_4-3\Preprocessing\170315_161220KI_4-3_b1_t150.tif';
o.TileFiles{2,1,1} = 'F:\In process\170315_161220KI_4-3\Preprocessing\170315_161220KI_4-3_b2_t150.tif';
o.TileFiles{3,1,1} = 'F:\In process\170315_161220KI_4-3\Preprocessing\170315_161220KI_4-3_b3_t150.tif';
o.TileFiles{4,1,1} = 'F:\In process\170315_161220KI_4-3\Preprocessing\170315_161220KI_4-3_b4_t150.tif';
o.TileFiles{5,1,1} = 'F:\In process\170315_161220KI_4-3\Preprocessing\170315_161220KI_4-3_b5_t150.tif';
o.TileFiles{6,1,1} = 'F:\In process\170315_161220KI_4-3\Preprocessing\170315_161220KI_4-3_SW_t150.tif';
o.CorrThresh=.4;    
o.ExtraCodes = {'Npy', 6, 3; 'Sst', 6, 5};
o.nExtraRounds=1;
o.OutputBaseName = '170315_161220KI_4-3';
o.SmoothSize = 1;


% make top-hat filtered tile files (for now a dummy function since this has
% already been done)
%o = iss_extract_and_filter(o);

% register tiles against their neighbors
o = iss_register(o);

% find spots and get their fluorescence values
[GlobalYX, SpotColors, Isolated] = iss_find_spots(o);
%GlobalYX = GoodGlobalYX; SpotColors = GoodSpotColors; Isolated = GoodIsolated;

% assign them to cells
[Genes, Codes, MaxScore, Intensity] = iss_call_spots(o, SpotColors, Isolated);

% this would make a figure before finishing 
% ScoreThresh = .85; ShowMe = (MaxScore>ScoreThresh);
% figure(100); 
% iss_make_figure(o, GlobalYX(ShowMe,:), Genes(ShowMe));

% find extra genes in final round (Sst and Npy)
[ExtraGlobalYX, ExtraGenes] = iss_single_genes(o);
%ExtraGlobalYX = ExtraGoodGlobalYX; ExtraGenes = ExtraGoodGene;

% now put everything together
FinalYX  = [GlobalYX ; ExtraGlobalYX];
FinalGenes = [Genes; ExtraGenes];
FinalMaxScore = [MaxScore ; ones(length(ExtraGenes),1)];

%% produce output figure
if o.Graphics
    figure(1)
    ShowMe = (FinalMaxScore>.9);
    iss_make_figure(o, FinalYX(ShowMe,:), FinalGenes(ShowMe));
   
end


%% now do DAPI processing

Dapi = imread(o.TileFiles{o.ReferenceRound,1,1}, o.DapiChannel);
[CellMap, CellYX] = iss_segment_dapi(o, Dapi);

figure(12); clf
ShowMe = (FinalMaxScore>.9);
iss_make_figure(o, FinalYX(ShowMe,:), FinalGenes(ShowMe));
% hold on; 
% plot(CellYX(:,2), CellYX(:,1), 'wo', 'markersize', 10);

%% call cells
load E:\PROOOJECTS\1_Neuron_mapping\SingleCell\1K\gSet
% if ~exist('gSet'); load gSet; end

SpotYX = FinalYX(FinalMaxScore>.9,:);
SpotGenes = FinalGenes(FinalMaxScore>.9);
SpotGenes(strcmp('Lphn2',SpotGenes))={'Adgrl2'}; % don't know why gene name changed

[pCellClass, pSpotCell] = iss_call_cells(o, SpotYX, SpotGenes, gSet, CellMap);


%% and plot pies
figure(2340976)
iss_pie_plot(o, CellYX, pCellClass);
%% save it
% savefig([o.OutputDirectory '\' o.OutputBaseName '.fig']);
% 
% % and save data
% save([o.OutputDirectory '\' o.OutputBaseName '.mat'], 'FinalYX', 'FinalGenes', 'FinalMaxScore', 'CellYX', 'CellMap');
