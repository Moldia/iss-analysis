% ISS image analysis workshop, 170614
% Xiaoyan
% tested on MATLAB R2016b

close all; clear;
cd SingleCell

%%
tilesize = 2000;
ntilesX = 2;
ntilesY = 1;
% 
% ndigits = 2;
tilepos = [0 0; 2000 0];

% mkdir('Stitched');
% mkdir('ParentCell');
% mkdir('CellBlobs');

%% relate blobs to parent cell (if not done in CellProfiler)
% Icell = cell(ntilesY,ntilesX);
% for i = 1:ntilesY
%     fprintf('%.2f percet read.\n', (i-1)/ntilesY*100);
%     for j = 1:ntilesX
%         t = ntilesX*(i-1)+j;
%         I = imread(['Segmentation\CellLabel_', paddigits(t,ndigits), '.tif']);
%         Icell(i,j) = {relabelimg(I)};
%     end  
% end  
% clear I
% fprintf('100 percet finished.\n');
% save ..\WS170614 Icell -append
load ..\WS170614 Icell

figure(1); clf;
set(gcf, 'name', 'Neighboring tiles', 'units', 'normalized', 'position', [.1 .2 .8 .6]);
imshow(label2rgb([Icell{1}, Icell{2}]));
title('labeled images of two neighboring tiles');
drawnow
hold on;

[~, pos] = getinsitudata('QT_0.41_0.004_details.csv');
pos = correctcoord(pos, 1);     % difference between zero indexing (Python) and non-zero indexing
cellid = zeros(length(pos),1);
for i = 1:ntilesY
    for j = 1:ntilesX
        in = readsinsqr(pos, tilesize*[j-1 i-1 j i]+[10000 4000 10000 4000]+.5);   % tile 46 starts at (10000, 4000)
        postmp = bsxfun(@minus, pos(in,:), tilesize*[j-1 i-1]+[10000 4000]);
        plot(postmp(:,1)+2000*(j-1), postmp(:,2), '+');
        postmp = readsinroi(postmp, Icell{i,j});
        cellid(in) = postmp;
    end
end

% writeblobwcell('QT_0.41_0.004_details.csv', cellid, 'QT_0.41_0.004');
% movefile('ParentCell\QT_0.41_0.004_details_wCell.csv', 'ParentCell\QT_0.41_0.004_details_wNewTileCell.csv');

%% relabel tile object images outputed from CP
% Ilabeled = cell(ntilesY,ntilesX);
% Ioutlines = cell(ntilesY,ntilesX);
% 
% for i = 1:ntilesY
%     fprintf('%.2f percet read.\n', (i-1)/ntilesY*100);
%     for j = 1:ntilesX
%         t = ntilesX*(i-1)+j;
%         I = imread(['Segmentation\NucleiLabel_', paddigits(t,ndigits), '.tif']);
%         Ivis = imread(['Outlines\overlayOutlines_' paddigits(t,ndigits) '.png']);
%         I = relabelimg(I);
%         Ilabeled(i,j) = {I};
%         Ioutlines(i,j) = {Ivis};
%     end
%     
% end  
% clear I Ivis
% fprintf('100 percet finished.\n');
% save ..\WS170614 Ilabeled Ioutlines -append
load ..\WS170614 Ilabeled Ioutlines

% stitched outline images for visualization
resizef = .5;
Ivis = resizestitch([ntilesX,ntilesY], tilesize, resizef, Ioutlines);

figure(1); clf;
set(gcf, 'name', 'Neighboring tiles', 'units', 'normalized', 'position', [.1 .2 .8 .6]);
imshow(Ivis);
title('segmented boundaries of two neighboring tiles');
drawnow

%% cell properties
cellprop = csvread('Cells.csv', 1);     % needs area, centroid position

% cell global position
if size(cellprop,2) > 3
    for i = 1:max(cellprop(:,1))
        try
            cellprop(cellprop(:,1)==i, end-1:end) = ...
                bsxfun(@plus, cellprop(cellprop(:,1)==i,end-1:end), tilepos(i,:));
        end
    end
else
    warning('Too few columns detected in cell property file.') 
end

% add metadata position column
cellprop = [cellprop(:,1:2), cellprop(:,1), cellprop(:,3:end)];

%% find objects cut by tiling lines
[CellCorrectionTable, CellProps, Ioutlines] = correct_edge_objects...
    (Ilabeled, Ioutlines, cellprop, tilepos, Ivis, resizef);
% save('Stitched\CellLookupTable.mat', 'CellCorrectionTable', 'CellProps');

% remake the resized outline image
Ivis = resizestitch([ntilesX,ntilesY], tilesize, resizef, Ioutlines);
% imwrite(Ivis, ['Stitched\Outlines_' num2str(resizef*100) '%.jpg']);
IvisHD = resizestitch([ntilesX,ntilesY], tilesize, 1, Ioutlines);
% imwrite(IvisHD, 'Stitched\Outlines_100%.jpg');

clear Ilabeled Ioutlines IvisHD

%% correct cell properties (area, centroid position)
%  properties should be from expanded objects ("real parents")
figure,imshow(Ivis); title('objects after correction')
hold on

% renumber cells
cellrenumber = renumberedgeobj...
    ([ntilesX,ntilesY], CellCorrectionTable, CellProps, 2);
[uCells, ~, idxCell] = unique(cellrenumber);
cCell = hist(idxCell, 1:length(uCells));
% fid = fopen('Stitched\UniqueCells.txt', 'w');
% fprintf(fid, '%d\n', uCells);
% fclose(fid);

CellPropsRenum = [uCells(cCell==1),...
    CellProps(ismember(idxCell, find(cCell==1)),3:6)];
plot(CellPropsRenum(:,end-1)*resizef, CellPropsRenum(:,end)*resizef, 'y.');

% merge cells and update properties
edgecells = uCells(cCell~=1);
for i = 1:length(edgecells)
    if edgecells(i) > 0
        edgecellProps = CellProps(cellrenumber==edgecells(i),3:6);
        % reassign tile number according to biggest part and sum up area
        [~, tilere] = max(edgecellProps(:,2));
        tile = edgecellProps(:,1);
        tilere = tile(tilere);
        atotal = sum(edgecellProps(:,2),1);
        plot(edgecellProps(:,3)*resizef,edgecellProps(:,4)*resizef,'-');
        % recalculate center of mass
        pos = edgecellProps(:,3:4).*repmat(edgecellProps(:,2)/atotal, 1, 2);
        try
            pos = sum(pos,1);
        catch
            pos = [nan, nan];
        end
        plot(pos(:,1)*resizef,pos(:,2)*resizef,'o');
        CellPropsRenum = [CellPropsRenum; [edgecells(i), tilere, atotal, pos]];
    end
end

[~,idxsort] = sort(CellPropsRenum(:,1));
CellPropsRenum = CellPropsRenum(idxsort,:);
% fid = fopen('Stitched\ExpandedCells.csv', 'w');
% fprintf(fid, 'CellID,metadata_position,area,global_x_pos,global_y_pos\n');
% fprintf(fid, '%d,%d,%d,%d,%d\n', reshape(CellPropsRenum',[],1));
% fclose(fid);
% save('Stitched\CellLookupTable.mat', 'CellPropsRenum', '-append');

%% blob cell relation
% parentcell = csvread('blobs.csv', 1);
parentcell = csvread('ParentCell\QT_0.41_0.004_details_wNewTileCell.csv', 1, 3);
columnCell = 6;

% original tile size different cell segmentation tile size
[~, pos] = getinsitudata('QT_0.41_0.004_details.csv');
pos = correctcoord(pos, 1);     % difference between zero indexing (Python) and non-zero indexing
tileid = ceil(pos/2000);
% original image 20 2000x2000 tiles in X direction
tileid = (tileid(:,2)-1)*20 + tileid(:,1);
tileid = tileid - 45;
tileid(tileid<0 | tileid>2) = 0;
parentcell(:,3) = tileid;

parentcellNew = renumberedgeobj...
    ([ntilesX,ntilesY], CellCorrectionTable, parentcell, columnCell);

% towrite = [(1:size(parentcell,1))', parentcell(:,3), parentcellNew]';
% fid = fopen('ParentCell\ParentCell_AllRefBlobs.csv', 'w');
% fprintf(fid,'blob_num, metadata_position, Parent_Cell\n');
% fprintf(fid, '%d,%d,%d\n', towrite(:));
% fclose(fid);                                                                    

clf; imshow(Ivis); hold on;
randcell = randi(max(parentcellNew), 100, 1);
gscatter((pos(ismember(parentcellNew, randcell),1)-10000)*resizef,...
    (pos(ismember(parentcellNew, randcell),2)-4000)*resizef,...
    parentcellNew(ismember(parentcellNew, randcell)));
legend off

%% label child blobs and get single cell profile
% % before QT
% parent = findparentcell('..\CP_170122_Seq\Decoding\beforeQT.mat',...
%     'blob_allbt', parentcellNew);
% csvwrite('ParentCell\ParentCell_beforeQT.csv', parent);
% writeblobwcell('..\CP_170122_Seq\Decoding\beforeQT_details.csv',...
%     parent, 'beforeQT');
% childblobs('ParentCell\beforeQT_details_wCell.csv',...
%     parent, CellPropsRenum, 'beforeQT');

% after QT
% parent = findparentcell('..\CP_170122_Seq\Decoding\QT_0.45_0.0001.mat',...
%     'blob_allqt', parentcellNew);
% csvwrite('ParentCell\ParentCell_QT_0.45.csv', parent);
writeblobwcell('QT_0.41_0.004_details.csv',...
    parentcellNew, 'QT_0.41');
childblobs('ParentCell\QT_0.41_details_wCell.csv',...
    parentcellNew, CellPropsRenum, 'QT_0.41');
