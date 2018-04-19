% cell calling on already cropped CA1
% modified from iss suite by Kenneth Harris
% Xiaoyan, 2018

%%
% roi cropped raw dapi image
cropped_cell_image = imread('C:\Users\markus.hilscher\Documents\MATLAB\schizo_models\ROI\170717_schizo_CA1Probes\1441-hippo\Ab_c1_ROI1.tif');
% name and position of spots in roi (coordinates match he above image)
[name, pos] = getinsitudata('C:\Users\markus.hilscher\Documents\MATLAB\schizo_models\ROI\170717_schizo_CA1Probes\1441-hippo\spots_ROI1.csv',...
    1);

% single cell data
load gSetCA1all.mat

%% segment DAPI
[CellMap, DapiBoundaries] = fnc_segment_dapi(cropped_cell_image);


%% probablistic cell calling and make pie charts
% change Lphn2 to Adgrl2
name(strcmp(name, 'Lphn2')) = {'Adgrl2'};

[CellYX ,pSpotCell ,pCellClass] =...
    fnc_call_cells_scClassPrior(CellMap, gSet, name, pos);

axis image

