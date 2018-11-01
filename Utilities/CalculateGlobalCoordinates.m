%% Calculate global coordinates based on the Python tiling file
%  for blob detection result
%  Xiaoyan, 2014-12-12


%% parameters
imageanalysis_file = 'input_example\blobs_ACTB.csv';
position_file = 'input_example\Tiled.csv';
transcript_name = 'ACTB';
output_filename = 'output_example\ACTB_decoded.csv';

%% calculate positions
tile_start_pos = getposition(position_file);

seqdata = csvread(imageanalysis_file,1);
tilenum = max(seqdata(:,3));
blobs = zeros(size(seqdata,1),4);

blobnum = 0;
for t = 1:tilenum
    tile_seq_data = seqdata(seqdata(:,3) == t,:);
    tile_blobnum = size(tile_seq_data,1);
    if tile_blobnum
        tile_x_start = mean(tile_start_pos(tile_start_pos(:,1) == t,2));
        tile_y_start = mean(tile_start_pos(tile_start_pos(:,1) == t,3));
        
        blobs(blobnum+1:blobnum+tile_blobnum,1:2) = ...
            bsxfun(@plus,tile_seq_data(:,5:6),[tile_x_start,tile_y_start]);
        blobs(blobnum+1:blobnum+tile_blobnum,4) = t;
        
        blobnum = blobnum+tile_blobnum;
    end
end
        
blobs_write = [repmat({'detection',transcript_name},size(seqdata,1),1),...
    num2cell(blobs)];
blobs_write = blobs_write';

%% plot
figure;
plot(blobs(:,1),blobs(:,2),'.');
set(gca,'YDir','reverse');
axis image

%% write
fid = fopen(output_filename,'w');
fprintf(fid,'letters,name,global_X_position,global_Y_position,parent_cell,tile_ID\n');
fprintf(fid,'%s,%s,%d,%d,%d,%d\n',blobs_write{:});
fclose(fid);