%deprecated 


load('..\CP_161031_Kv21\Stitched\CellLookupTable.mat');
tilepos = getcsvposition('..\ZENout\Tiled_161031_Kv21.csv');

tile_size = tilepos(2,2)-tilepos(1,2);
tiles_x = range(tilepos(:,2))/tile_size + 1;
tiles_y = range(tilepos(:,3))/tile_size + 1;

digit_num = 3;

mkdir('Stitched\');
mkdir('ParentCell\')

% correct blob celll relation
parentcell = csvread('Sst.csv', 1);
column_cell = 10;

for i = 1:tiles_y
    for j = 1:tiles_x        
        tilenum = tiles_x*(i-1)+j;
        from = parentcell(parentcell(:,3)==tilenum,column_cell);

        if ~isempty(from)
            to = CellTable{i,j};
            for k = 1:size(from,1)
                if from(k)
                    try
                        from(k) = to(to(:,1)==from(k),2);
                    catch
                        from(k) = -1;   % cell not found in the object image;
                    end
                end
            end
            parentcell(parentcell(:,3)==tilenum,column_cell) = from;
        end
    end
end

towrite = [(1:size(parentcell,1))', parentcell(:,3), parentcell(:,column_cell)]';
fid = fopen('ParentCell\ParentCell_AllRefBlobs.csv', 'w');
fprintf(fid,'blob_num, metadata_position, Parent_Cell\n');
fprintf(fid, '%d,%d,%d\n', towrite(:));
fclose(fid);

%% label child blobs
%  before QT

load('Decoding\beforeQT.mat', 'blob_allbt');
parentcell_bt = parentcell(blob_allbt,column_cell);
csvwrite('ParentCell\ParentCell_beforeQT.csv', parentcell_bt);

blob = importdata('Decoding\beforeQT_details.csv');
fid = fopen('ParentCell\beforeQT_details_wCell.csv','w');
formspec = [];
for i = 1:size(blob.textdata,2)+1
    formspec = [formspec, '%s,'];
end
formspec = [formspec(1:end-1), '\n'];
header = [blob.textdata(1,:),{'Parent_Cell'}];
fprintf(fid,formspec,header{:});
blobwrite = [blob.textdata(2:end,1:2),num2cell([blob.data,double(parentcell_bt)])]';
formspec = '%s,%s,';
for i = 1:size(blob.textdata,2)-1
    formspec = [formspec, '%d,'];
end
formspec = [formspec(1:end-1), '\n'];
fprintf(fid,formspec,blobwrite{:});
fclose(fid);

%% label child blobs
%  after QT

load('Decoding\QT_0.5_0.mat', 'blob_allqt');
parentcell_qt = parentcell(blob_allqt,column_cell);
csvwrite('ParentCell\ParentCell_QT_0.5.csv', parentcell_qt);

blob = importdata('Decoding\QT_0.5_0_details.csv');
fid = fopen('ParentCell\QT_0.5_details_wCell.csv','w');
formspec = [];
for i = 1:size(blob.textdata,2)+1
    formspec = [formspec, '%s,'];
end
formspec = [formspec(1:end-1), '\n'];
header = [blob.textdata(1,:),{'Parent_Cell'}];
fprintf(fid,formspec,header{:});
blobwrite = [blob.textdata(2:end,1:2),num2cell([blob.data,double(parentcell_qt)])]';
formspec = '%s,%s,';
for i = 1:size(blob.textdata,2)-1
    formspec = [formspec, '%d,'];
end
formspec = [formspec(1:end-1), '\n'];
fprintf(fid,formspec,blobwrite{:});
fclose(fid);
