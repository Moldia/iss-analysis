%% convert blob-cell relation to cellular profile
%  needs both decoding and parent cell file
%  Xiaoyan, last update 161106
% deprecated

files = {'QT_0.45_1e-05','beforeQT'};
filenames = {'QT_0.45','beforeQT'};
mkdir('CellBlobs');

for m = 1:length(files)
    
    blobfile = ['..\CP_161014\Decoding\' files{m} '_details.csv'];
    [names, pos] = getinsitudata_f(blobfile);
    [name_uni, ~, idx_name] = unique(names);
    
    cells = csvread(['ParentCell\ParentCell_' filenames{m} '.csv']);
    [cell_uni, ~, idx_cell] = unique(cells);
    
    Count = zeros(length(cell_uni),length(name_uni));
    
    for i = 1:length(name_uni)
        temp = idx_cell(idx_name==i);
        Count(:,i) = hist(temp,1:length(cell_uni));
    end
    
    %% cell properties
    cellprops = csvread('Stitched\ExpCells.csv', 1);
    
    cellnum = zeros(length(cell_uni),1);
    for i = 1:length(cell_uni)
        if cell_uni(i) > 0
            cellnum(i) = find(cellprops(:,1)==cell_uni(i));
        else
            cellnum(i) = 0;
        end
    end
    
    if  cell_uni(1) == -1
        if ~cell_uni(2)
            cells_wblobs = [nan(2, 4); cellprops(cellnum(3:end),2:end)];
        else
            cells_wblobs = [nan(1, 4); cellprops(cellnum(2:end),2:end)];
        end
    else
        if ~cell_uni(1)
            cells_wblobs = [nan(1, 4); cellprops(cellnum(2:end),2:end)];
        else
            cells_wblobs = cellprops(cellnum,2:end);
        end
    end
    cells_wblobs = [cell_uni, cells_wblobs];
    
    fid = fopen(['Stitched\ExpCells_wBlobs_' filenames{m} '.csv'], 'w');
    fprintf(fid,'CellID, metadata_position,area,global_x_pos,global_y_pos\n');
    fprintf(fid,'%d,%d,%d,%d,%d\n',reshape(cells_wblobs',[],1));
    fclose(fid);
    
    %% cell and child blobs, raw, including NNNN and cell id = 0
    fid = fopen(['CellBlobs\CellBlobs_raw_' filenames{m} '.csv'],'w');
    
    name_write = [{'CellID'}, name_uni'];
    formspec = '%s';
    for i = 1:length(name_uni)
        formspec = [formspec, ',%s'];
    end
    fprintf(fid,formspec,name_write{:});
    fprintf(fid,',centroid_x,centroid_y\n');
    
    mwrite = [cell_uni,Count,cells_wblobs(:,end-1:end)];
    formspec = [];
    for i = 1:length(name_uni)+1
        formspec = [formspec, '%d,'];
    end
    formspec = [formspec, '%d,%d\n'];
    fprintf(fid,formspec,mwrite');
    
    fclose(fid);
    
    %% cell and child blobs, remove NNNN and cell id = 0 and cell id = -1
    fid = fopen(['CellBlobs\CellBlobs_' filenames{m} '.csv'],'w');
    
    idx_NNNN = strcmp(name_uni,'NNNN');
    rowkeep = sum(Count(:,~idx_NNNN),2) ~= 0;
    if cell_uni(1) == -1
        rowkeep(1) = false;
        if ~cell_uni(2)
            rowkeep(2) = false;
        end
    elseif ~cell_uni(1)
        rowkeep(1) = false;
    end
    
    name_write = [{'cellID'}, name_uni(~idx_NNNN)'];
    formspec = '%s';
    for i = 1:length(name_uni)-1
        formspec = [formspec, ',%s'];
    end
    fprintf(fid,formspec,name_write{:});
    fprintf(fid,',centroid_x,centroid_y\n');
    
    mwrite = [cell_uni(rowkeep),Count(rowkeep,~idx_NNNN),cells_wblobs(rowkeep,end-1:end)];
    formspec = [];
    for i = 1:length(name_uni)
        formspec = [formspec, '%d,'];
    end
    formspec = [formspec, '%d,%d\n'];
    fprintf(fid,formspec,mwrite');
    
    fclose(fid);
end


%% visualization
I = imread('Stitched\Outlines_20%.png');
figure
imshow(I)
hold on
cells = csvread('CellBlobs\CellBlobs_beforeQT.csv',1);
plot(cells(:,end-1)*.2,cells(:,end)*.2,'yp')
cells = csvread('CellBlobs\CellBlobs_QT_0.45.csv',1);
plot(cells(:,end-1)*.2,cells(:,end)*.2,'wp')