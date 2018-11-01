% read position from .czi metadata and write CP input file
% only for tile series images
% Xiaoyan, 2017

%% input
czifile = 'F:\CARTANA\spinal code\2018-05-15-SC_b1_12\180515_mSC_b1.czi';

tif_file_basename_prefix = 'base';
tif_file_scenename_prefix = '_s';
tif_file_channelname_prefix = 'c';
tif_file_suffix = '_ORG.tif';

channel_name = {'DAPI', 'anchor', 'A', 'G', 'C' 'T'};

outputfolder = 'G:\Lotta\20180602\ZENout';

%%
% get tile position from metadata
[nSeries, nSerieswPos, ~, ~, xypos] = get_ome_tilepos(czifile);

% exported files
imdirs = cell(nSerieswPos, 4);   % default: 4 bases, 6 channels
imfiles = cell(nSerieswPos, 4, 6);
for s = 1:nSerieswPos
    for b = 1:4
        imdirs{s,b} = [tif_file_basename_prefix, num2str(b)];
        for c = 1:6
            imfiles{s,b,c} = [tif_file_scenename_prefix, paddigits(s,2),...
                tif_file_channelname_prefix, num2str(c), tif_file_suffix];
        end
    end
end

% =========================================================================
% JUST BECAUSE NOT ALL BASES HAVE THE SAME FILE NAME STRUCTURE
imfiles(:,1,:) = cellfun(@(v)...
    strrep(v, tif_file_scenename_prefix, 'new-Orthogonal Projection-01_s'),...
    imfiles(:,1,:), 'uni', 0);
% =========================================================================

% per hyb data
towrite = [];
header = channel_name;
metadata = num2cell([reshape(1:45, [], 1), xypos]);
for b = 1:4
    perbase = [];
    for c = 1:6
        perbase = [perbase, imdirs(:,b), imfiles(:,b,c)];
    end
    towrite = [towrite; [metadata, num2cell(repmat(b, nSerieswPos, 1)), perbase]];
end

% add reference channels (default 1st (DAPI) and 2nd (Cy7) in hyb1)
refchannels = 1:2;
for i = 1:numel(refchannels)
    towrite = [towrite,...
        repmat([imdirs(:,1), imfiles(:,1,i)], 4, 1)];
    header = [header, ['General_' channel_name{refchannels(i)}]];
end
[~, idx] = sort(cell2mat(towrite(:,1)));
towrite = towrite(idx,:)';

% file header
header = [strcat('Image_PathName_', header); strcat('Image_FileName_', header)];
header = [{'Metadata_position','Tile_xPos','Tile_yPos','Hyb_step'}, reshape(header, 1, [])];

% write file
mkdir(outputfolder)
fid = fopen(fullfile(outputfolder, 'CPinput.csv'), 'w');
fprintf(fid, lineformat('%s', length(header)), header{:});
fprintf(fid, ['%d,%d,%d,%d,' lineformat('%s', length(header)-4)], towrite{:});
fclose(fid);
