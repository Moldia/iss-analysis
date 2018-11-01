% re-relate blob and already edge-corrected parent cells
% Xiaoyan, 2017

tilepos = getcsvtilepos('..\Preprocessing\Stitched\Tiled_170206_Cell.csv');
tilesize = tilepos(2,2)-tilepos(1,2);
ntilesX = range(tilepos(:,2))/tilesize + 1;
ntilesY = range(tilepos(:,3))/tilesize + 1;

mkdir('ParentCell')

load('..\CP_170206_Cell\Stitched\CellLookupTable.mat');

SW = {'Sst', 'Npy'};

for i = 1:length(SW)
    % correct blob cell relation
    parentcell = csvread([SW{i}, '.csv'], 1);
    columnCell = 7;
    
    parentcellNew = renumberedgeobj...
        ([ntilesX,ntilesY], CellCorrectionTable, parentcell, columnCell);
    
    towrite = [(1:size(parentcell,1))', parentcell(:,3), parentcellNew]';
    fid = fopen(['ParentCell\ParentCell_', SW{i}, '.csv'], 'w');
    fprintf(fid,'blob_num, metadata_position, Parent_Cell\n');
    fprintf(fid, '%d,%d,%d\n', towrite(:));
    fclose(fid);
    
    % label child blobs and get single cell profile
    parent = findparentcell(['Decoding_', SW{i}, '\beforeQT.mat'],...
        'blob_allbt', parentcellNew);
    csvwrite(['ParentCell\ParentCell_', SW{i}, '.csv'], parent);
    writeblobwcell(['Decoding_', SW{i}, '\beforeQT_details.csv'],...
        parent, SW{i});
end

