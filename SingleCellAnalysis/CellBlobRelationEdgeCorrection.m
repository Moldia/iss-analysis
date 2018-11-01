% correct cell objects cut by tiling line
% re-relate blob and parent
% Xiaoyan, 2017


% parpool
tilepos = getcsvtilepos('..\Tiled_170526.csv');

tilesize = tilepos(2,2)-tilepos(1,2);
ntilesX = range(tilepos(:,2))/tilesize + 1;
ntilesY = range(tilepos(:,3))/tilesize + 1;

ndigits = ceil(log10(ntilesX*ntilesY));

mkdir('Stitched');
mkdir('ParentCell')
mkdir('CellBlobs');

%% relabel tile object images outputed from CP
Ilabeled = cell(ntilesY,ntilesX);
Ioutlines = cell(ntilesY,ntilesX);

for i = 1:ntilesY
    fprintf('%.2f percet read.\n', (i-1)/ntilesY*100);
    parfor j = 1:ntilesX
        t = ntilesX*(i-1)+j;
        I = imread(['CP_170526_Cell\Segmentation\NucleiLabel_', paddigits(t,ndigits), '.tif']);
        Ivis = imread(['CP_170526_Cell\Outlines\overlayOutlines_' paddigits(t,ndigits) '.png']);
        I = relabelimg(I);
        Ilabeled(i,j) = {I};
        Ioutlines(i,j) = {Ivis};
    end
    
end  
clear I Ivis
fprintf('100 percet finished.\n');

% stitched outline images for visualization
resizef = .5;
Ivis = resizestitch([ntilesX,ntilesY], tilesize, resizef, Ioutlines);

%% relate blobs to parent cell (if not done in CellProfiler)
Icell = cell(ntilesY,ntilesX);
for i = 1:ntilesY
    fprintf('%.2f percet read.\n', (i-1)/ntilesY*100);
    parfor j = 1:ntilesX
        t = ntilesX*(i-1)+j;
        I = imread(['CP_170526_Cell\Segmentation\CellLabel_', paddigits(t,ndigits), '.tif']);
        Icell(i,j) = {relabelimg(I)};
    end  
end  
clear I
fprintf('100 percet finished.\n');

[~, pos] = getinsitudata('..\QT_0.4_0.004_details.csv');
pos = correctcoord(pos, 1);     % difference between zero indexing (Python) and non-zero indexing
cellid = zeros(length(pos),1);
for i = 1:ntilesY
    for j = 1:ntilesX
        in = readsinsqr(pos, tilesize*[j-1 i-1 j i]+.5);
        postmp = bsxfun(@minus, pos(in,:), tilesize*[j-1 i-1]);
        postmp = readsinroi(postmp, Icell{i,j});
        cellid(in) = postmp;
    end
end
writeblobwcell('..\QT_0.4_0.004_details.csv', cellid, 'ParentCell\QT_0.4_0.004');

%% cell properties
cellprop = csvread('CP_170526_Cell\Cells.csv', 1);     % needs area, centroid position

% cell global position
if size(cellprop,2) > 3
    for i = 1:max(cellprop(:,3))
        try
            cellprop(cellprop(:,3)==i, end-1:end) = ...
                bsxfun(@plus, cellprop(cellprop(:,3)==i,end-1:end), tilepos(i,2:3));
        end
    end
else
    warning('Too few columns detected in cell property file.') 
end

%% find objects cut by tiling lines
[CellCorrectionTable, CellProps, Ioutlines] = correct_edge_objects...
    (Ilabeled, Ioutlines, cellprop, tilepos, Ivis, resizef);
save('Stitched\CellLookupTable.mat', 'CellCorrectionTable', 'CellProps');

% remake the resized outline image
Ivis = resizestitch([ntilesX,ntilesY], tilesize, resizef, Ioutlines);
imwrite(Ivis, ['Stitched\Outlines_' num2str(resizef*100) '%.jpg']);
IvisHD = resizestitch([ntilesX,ntilesY], tilesize, 1, Ioutlines);
imwrite(IvisHD, 'Stitched\Outlines_100%.jpg');

clear Ilabeled Ioutlines IvisHD

%% correct cell properties (area, centroid position)
%  properties should be from expanded objects ("real parents")
figure,imshow(Ivis)
hold on

% renumber cells
cellrenumber = renumberedgeobj...
    ([ntilesX,ntilesY], CellCorrectionTable, CellProps, 2);
[uniCell, ~, idxCell] = unique(cellrenumber);
countCell = hist(idxCell, 1:length(uniCell));
fid = fopen('Stitched\UniqueCells.txt', 'w');
fprintf(fid, '%d\n', uniCell);
fclose(fid);

CellPropsRenum = [uniCell(countCell==1),...
    CellProps(ismember(idxCell, find(countCell==1)),3:6)];
plot(CellPropsRenum(:,end-1)*resizef, CellPropsRenum(:,end)*resizef, '.');

% merge cells and update properties
edgecells = uniCell(countCell~=1);
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
fid = fopen('Stitched\ExpandedCells.csv', 'w');
fprintf(fid, 'CellID,metadata_position,area,global_x_pos,global_y_pos\n');
fprintf(fid, '%d,%d,%d,%d,%d\n', reshape(CellPropsRenum',[],1));
fclose(fid);
save('Stitched\CellLookupTable.mat', 'CellPropsRenum', '-append');

%% blob cell relation
% parentcell = csvread('blobs.csv', 1);
parentcell = csvread('QT_0.4_0.004_details_wCell.csv', 1, 3);
columnCell = 6;

% original tile size different cell segmentation tile size
tileid = ceil(pos/2000);
tileid = (tileid(:,2)-1)*ntilesX + tileid(:,1);
parentcell(:,3) = tileid;

parentcellNew = renumberedgeobj...
    ([ntilesX,ntilesY], CellCorrectionTable, parentcell, columnCell);

towrite = [(1:size(parentcell,1))', parentcell(:,3), parentcellNew]';
fid = fopen('ParentCell\ParentCell_AllRefBlobs.csv', 'w');
fprintf(fid,'blob_num, metadata_position, Parent_Cell\n');
fprintf(fid, '%d,%d,%d\n', towrite(:));
fclose(fid);                                                                    

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
writeblobwcell('..\QT_0.4_0.004_details.csv',...
    parentcellNew, 'QT_0.4');
childblobs('ParentCell\QT_0.4_details_wCell.csv',...
    parentcellNew, CellPropsRenum, 'QT_0.4');

