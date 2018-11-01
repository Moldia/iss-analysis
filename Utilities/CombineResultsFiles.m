%% combine sequencing and detection files
%  and CellBlobs file
%  Xiaoyan, 2016-9-9

%% sequencing/ detection results
files = {'CP_170206_Cell\ParentCell\QT_0.45_details_wCell.csv';
    'CP_170208_SWCell\ParentCell\Npy_details_wCell.csv';
    'CP_170208_SWCell\ParentCell\Sst_details_wCell.csv'};

type = {'sequencing', 'sequencing', 'sequencing'};

towrite = [];
if length(unique(type))==1
    for i = 1:length(files)
        filedata = importdata(files{i});
        header = filedata.textdata(1,:);
        towrite = [towrite; [filedata.textdata(2:end,1:2), num2cell(filedata.data)]];
    end
    towrite = towrite';
    fid = fopen('CombinedResults_wCell.csv', 'w');
    fprintf(fid, lineformat('%s', length(header)), header{:});
    fprintf(fid,['%s,%s,', lineformat('%d', length(header)-2)], towrite{:});
    fclose(fid);    
else
    for i = 1:length(files)
        fid = fopen(files{i});
        if strcmp(type{i},'sequencing')
            filedata = textscan(fid, '%s%s%f%f%d%f%f%*[^\n]','headerlines',1,'delimiter',',');
            temp = cellfun(@(v) num2cell(v),filedata(3:7),'uni',0);
            towrite = [towrite; [[filedata{1:2}],[temp{:}]]];
        else
            filedata = textscan(fid, '%s%s%f%f%d%f%*[^\n]','headerlines',1,'delimiter',',');
            temp = cellfun(@(v) num2cell(v),filedata([3,4,6]),'uni',0);
            towrite = [towrite; [[filedata{1:2}],[temp{:}],num2cell(ones(length(filedata{1}),2))]];
            
        end
        fclose(fid);
    end
    towrite = towrite';
    fid = fopen('CombinedResults.csv', 'w');
    fprintf(fid,'letters,name,global_X_pos,global_Y_pos,tile_ID,general_stain_min,seq_quality_min\n');
    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d\n', towrite{:});
    fclose(fid);
end

%% CellBlobs files
mkdir('CellBlobs');
mkdir('Stitched');
load('CP_170206_Cell\Stitched\CellLookupTable.mat');
parent = importdata('CombinedResults_wCell.csv');
parent = parent.data(:,end);
childblobs('CombinedResults_wCell.csv',...
    parent, CellPropsRenum, 'Combined');
