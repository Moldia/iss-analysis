czi_file = '';  % full path
output_prefix = ''; % same as specified in ZEN while exporting original tile images
output_folder = ''; % where the exported images are

%%
[nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
    get_ome_tilepos(czi_file);

% get only file name
[~, output_prefix] = fileparts(output_prefix);

% tile position configuration file
csvwrite(fullfile(outdir, ['tile_coordinates_' output_prefix '.csv']), xypos);

% assumes ZEN output format
for c = 1:nChannels


    fid = fopen(fullfile(outdir, ['TileConfiguration_', output_prefix, '.txt']), 'w');
    fprintf(fid,'# Define the number of dimensions we are working on\n');
    fprintf(fid,'dim = 2\n\n# Define the image coordinates\n');
    for t = 1:nSerieswPos
        fprintf(fid,...
            '%s_c%dm%s_ORG.tif; ; (%.1f, %.1f)\n',...
            fullfile(output_folder, output_prefix),...
            c, paddigits(t, floor(c/10)), xypos(t,1), xypos(t,2));
    end
    fclose(fid);
end

