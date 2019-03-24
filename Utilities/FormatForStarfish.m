image_directory = 'F:\2nd_40x_week6.5\181208-original-MIP-Stitched\MIP_stitch_lessStack\Miped\tissue1';
image_files = {'base1(1).czi'
    'base2(1).czi'
    'base3(1).czi'
    'base4(1).czi'
    'base5(1).czi'};
output_directory = 'Preprocessing';
nChannels = 6;
reference_base = 1;
reference_channel = 1;

% [~, nFiles] = get2Dtiles(image_directory, image_files, output_directory);
% align2Dtiles(fullfile(output_directory, '2DTiles'), image_files, output_directory, nChannels, max(nFiles), 1, 1);
% stitch2DtilesISS(fullfile(output_directory, 'Aligned2DTiles'), image_files, output_directory, nChannels, 1, [.3 .6], 0);
% stitch2DtilesMIST(fullfile(output_directory, 'Aligned2DTiles'), image_files, output_directory, nChannels, 1, 0);
stitch2DtilesMIST(fullfile(output_directory, 'Aligned2DTiles'), image_files, output_directory, nChannels, 2, 1);

%% get MIPs or make EDFs from z stacks
function [outdir, nFiles] = get2Dtiles(indir, files, outdir)
outdir = fullfile(outdir, '2DTiles');
mkdir(outdir);

nFiles = zeros(numel(files));

for s = 1:numel(files)     % series
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    
    % Construct a Bio-Formats reader decorated with the Memoizer wrapper
    bfreader = loci.formats.Memoizer(bfGetReader(), 0);
    % Initialize the reader with an input file
    bfreader.setId(fullfile(imdir, imfile));
    
    [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
        get_ome_tilepos(bfreader);
    scene = nSeries/nSerieswPos;
    
    bfreader.close();
    
    nFiles(s) = nSerieswPos;
    
    for c = 1:nChannels
        
        parfor t = 1:nSerieswPos
            fprintf('Channel %d. Tile %d / %d\n',...
                c, t, nSerieswPos);
            
            % Initialize a new reader per worker as Bio-Formats is not thread safe
            bfreader2 = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % Initialization should use the memo file cached before
            bfreader2.setId(fullfile(indir, imfile));
            
            bfreader2.setSeries(scene*t-1);
            
            if nZstacks == 1
                % only one z plane in the file
                iPlane = bfreader2.getIndex(0, c-1, 0)+1;
                I = bfGetPlane(bfreader2, iPlane);
                
                imwrite(I, fullfile(outdir,...
                    [output_prefix, '_t', num2str(t) '.tif']),...
                    'tiff', 'writemode', 'append');
            else
                I = cell(nZstacks,1);
                for z = 1:nZstacks
                    iPlane = bfreader2.getIndex(z-1, c-1, 0)+1;
                    I{z} = bfGetPlane(bfreader2, iPlane);
                end
                % focus stacking
                IFS = fstack_modified(I);
                imwrite(IFS, fullfile(outdir,...
                    [output_prefix, '_t', num2str(t) '.tif']),...
                    'tiff', 'writemode', 'append');
            end
            bfreader2.close();
            
        end
        datetime
    end
    
    % tile position configuration file
    csvwrite(fullfile(outdir, ['tile_coordinates_' output_prefix '.csv']), xypos);
    
    fid = fopen(fullfile(outdir, ['TileConfiguration_', output_prefix, '.txt']), 'w');
    fprintf(fid,'# Define the number of dimensions we are working on\n');
    fprintf(fid,'dim = 2\n\n# Define the image coordinates\n');
    for t = 1:nSerieswPos
        fprintf(fid,...
            '%s_t%d.tif; ; (%.1f, %.1f)\n',...
            output_prefix, t, xypos(t,1), xypos(t,2));
    end
    fclose(fid);
end
end

%% align all tiles to the reference round
function align2Dtiles(indir, files, outdir, nChannels, nSeries, ref_round, ref_channel)
outdir = fullfile(outdir, 'Aligned2DTiles');
mkdir(outdir);

parfor t = 1:nSeries
    t
    
    for s = [ref_round, setdiff(1:numel(files), ref_round)]
        imfile = files{s};
        [~, output_prefix] = fileparts(imfile);
        
        if s == ref_round
            for c = 1:nChannels
                I = imread(fullfile(indir, [output_prefix '_t' num2str(t) '.tif']), c);
                imwrite(I, fullfile(outdir, [output_prefix '_t' num2str(t) '.tif']), 'WriteMode', 'append');
                if c == ref_channel
                    refimg = I;
                end
            end
        else
            for c = [ref_channel, setdiff(1:nChannels, ref_channel)]
                I = imread(fullfile(indir, [output_prefix '_t' num2str(t) '.tif']), c);
                if c == ref_channel
                    floatimg = I;
                    
                    % 0.5px-resolution registration
                    [tform, Greg] = dftregistration(fft2(refimg), fft2(floatimg), 2);
                end
                
                % but pixel-resolution transformation (no rotaiton)
                I = padimg(I, round(tform(4)), round(tform(3)), 'NW');
                I = padimg(I, 2048-size(I,2), 2048-size(I,1));
                imwrite(I, fullfile(outdir, [output_prefix '_t' num2str(t) '.tif']), 'WriteMode', 'append');
            end
        end
    end
end
end

%% iss_suite stitching, modified from @iss/register
function stitch2DtilesISS(indir, files, outdir, nChannels, ref_channel, corr_threshold, ref_file)

if nargin < 6
    corr_threshold = [.3 .6];
end

outdir_top = outdir;
if nargin < 7 || ref_file==0
    ref_file = 0;
    outdir = fullfile(outdir, 'Stitched2DTiles_ISS_Individual');
else
    outdir = fullfile(outdir, ['Stitched2DTiles_ISS_Ref' num2str(ref_file)]);
end

mkdir(outdir);

s_order = 1:numel(files);
if ref_file
    s_order = [ref_file, setdiff(s_order, ref_file)];
end

for s = s_order
    
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    
    xypos = csvread(fullfile(outdir_top, '2DTiles', ['tile_coordinates_' output_prefix '.csv']));
    
    
    [gridorder, tilefiles] = tilepos2grid(xypos,...
        [output_prefix '_t'], '.tif', 0);
    
    [nY, nX] = size(tilefiles);
    nTiles = nY*nX;
    NonemptyTiles = find(~strcmp(tilefiles, 'empty.tif'))';
    EmptyTiles = strcmp(tilefiles, 'empty.tif');
    
    if ~ref_file || s == ref_file
        
        try
            TileOrigin = csvread(fullfile(outdir, ['tile_origin_' output_prefix, '.csv']));
        catch
            
            RefImages = zeros(2048, 2048, nY, nX, 'uint16');
            
            for t = NonemptyTiles(:)'
                [y,x] = ind2sub([nY nX], t);
                if mod(t,10)==0; fprintf('Loading tile %d reference image\n', t); end
                Im = imread(fullfile(indir, tilefiles{y,x}), ref_channel);
                
                RefImages(:,:,t) = imfilter(Im, fspecial('disk', 1));
            end
            
            VerticalPairs = zeros(0,2);
            HorizontalPairs = zeros(0,2);
            vShifts = zeros(0,2);
            hShifts = zeros(0,2);
            ccv = zeros(0,1);
            cch = zeros(0,1);
            
            for t = NonemptyTiles
                [y,x] = ind2sub([nY nX], t);
                
                % align to south neighbor
                if y<nY && ~EmptyTiles(t+1)
                    [shift, cc, RefFFTStore] = ImRegFft2(RefImages(:,:,t), RefImages(:,:,t+1), corr_threshold, 200^2);
                    if all(isfinite(shift))
                        VerticalPairs = [VerticalPairs; t, t+1];
                        vShifts = [vShifts; shift];
                        ccv = [ccv; cc];
                    end
                    fprintf('%d, %d, down: shift %d %d, cc %f\n', y, x, shift, cc);
                end
                
                % align to east neighbor
                if x<nX && ~EmptyTiles(t+nY)
                    if isvarname(RefFFTStore)
                        [shift, cc] = ImRegFft2(RefImages(:,:,t), RefImages(:,:,t+nY), corr_threshold, 200^2);
                    else
                        [shift, cc] = ImRegFft2(RefImages(:,:,t), RefImages(:,:,t+nY), corr_threshold, 200^2);
                    end
                    if all(isfinite(shift))
                        HorizontalPairs = [HorizontalPairs; t, t+nY];
                        hShifts = [hShifts; shift];
                        cch = [cch; cc];
                    end
                    fprintf('%d, %d, right: shift %d %d, cc %f\n', y, x, shift, cc);
                    
                end
            end
            
            
            % solve a set of linear equations for each shift
            M = zeros(nTiles, nTiles);
            c = zeros(nTiles, 2);
            for i=1:size(VerticalPairs,1)
                if isnan(vShifts(i,1)); continue; end
                t1 = VerticalPairs(i,1);
                t2 = VerticalPairs(i,2);
                M(t1,t1) = M(t1,t1)+1;
                M(t1,t2) = M(t1,t2)-1;
                c(t1,:) = c(t1,:) - vShifts(i,:);
                M(t2,t2) = M(t2,t2)+1;
                M(t2,t1) = M(t2,t1)-1;
                c(t2,:) = c(t2,:) + vShifts(i,:);
            end
            
            for i=1:size(HorizontalPairs,1)
                if isnan(hShifts(i,1)); continue; end
                t1 = HorizontalPairs(i,1);
                t2 = HorizontalPairs(i,2);
                M(t1,t1) = M(t1,t1)+1;
                M(t1,t2) = M(t1,t2)-1;
                c(t1,:) = c(t1,:) - hShifts(i,:);
                M(t2,t2) = M(t2,t2)+1;
                M(t2,t1) = M(t2,t1)-1;
                c(t2,:) = c(t2,:) + hShifts(i,:);
            end
            
            % anchor one of the tiles to a fixed coordinate
            Huge = 1e6;
            TileDistFromCenter = abs(mod(0:nTiles-1, nY)-nY/2) + ...
                abs(floor((0:nTiles-1)/nY)-nX/2);
            [~, HomeTile] = min(TileDistFromCenter(:)./~EmptyTiles(:));
            %sub2ind([nY nX], ceil(nY/2), ceil(nX/2));
            M(nTiles+1,HomeTile) = 1;
            c(nTiles+1,:) = [Huge, Huge];
            
            Tiny = 1e-4; % for regularization
            TileOffset0 = (M+Tiny*eye(nTiles+1, nTiles))\c;
            
            % find tiles that are connected to the home tile
            AlignedOK = (TileOffset0(:,1)>Huge/2);
            TileOffset1 = nan(nTiles, 2);
            TileOffset1(AlignedOK,:) = TileOffset0(AlignedOK,:)-Huge;
            
            % RefPos(t,1:2) is origin of reference tile
            RefPos = bsxfun(@minus, TileOffset1, nanmin(TileOffset1))+1;
            
            % tile origin(t,1:2,r)
            TileOrigin =  round(RefPos);
            csvwrite(fullfile(outdir, ['tile_origin_' output_prefix, '.csv']), TileOrigin);
            
            clear RefImages
        end
    end
    
    % stitch image and write
    for c = 1:nChannels
        MaxTileLoc = max(TileOrigin);
        BigIm = zeros(ceil((MaxTileLoc + 2048)), 'uint16');
        for t = NonemptyTiles
            MyOrigin = TileOrigin(t,:);
            if ~isfinite(MyOrigin(1)); continue; end
            if mod(t,10)==0; fprintf('Loading channel %d tile %d image\n', c, t); end
            LocalIm = imread(fullfile(indir, tilefiles{t}), c);
            BigIm(floor(MyOrigin(1))+(1:2048), floor(MyOrigin(2))+[1:2048]) ...
                = imresize(LocalIm, 1);
        end
        imwrite(BigIm, fullfile(outdir, [output_prefix, '_c', num2str(c) '.tif']));
    end
    
end
end

%% MIST stitching
function stitch2DtilesMIST(indir, files, outdir, nChannels, ref_channel, ref_file)

outdir_top = outdir;
if nargin < 6 || ref_file==0
    ref_file = 0;
    outdir = fullfile(outdir, 'Stitched2DTiles_MIST_Individual');
else
    outdir = fullfile(outdir, ['Stitched2DTiles_MIST_Ref' num2str(ref_file)]);
end

mkdir(outdir);

s_order = 1:numel(files);
if ref_file
    s_order = [ref_file, setdiff(s_order, ref_file)];
end

save_stitched_image = 1;

imwrite(zeros(2048, 'uint16'), [indir '\empty.tif']);

for s = s_order
    
    imfile = files{s};
    [~, output_prefix] = fileparts(imfile);
    
    xypos = csvread(fullfile(outdir_top, '2DTiles', ['tile_coordinates_' output_prefix '.csv']));
    
    [gridorder, tilefiles] = tilepos2grid(xypos,...
        [output_prefix '_t'], '.tif', 0);
    
    output_prefix = [output_prefix '_'];

    if ref_file && s~=ref_file
        % use exisiting tramsformation
        imfile = files{ref_file};
        [~, output_prefix_ref] = fileparts(imfile);
        output_prefix_ref = [output_prefix_ref '_'];
        
        assemble_from_metadata = 1;     % apply the same to other images
        ref_transform = fullfile(outdir, [output_prefix_ref 'metadata-' num2str(ref_channel) '.mat']);
        
        for c = [ref_channel, setdiff(1:nChannels, ref_channel)]            
            copyfile(ref_transform, fullfile(outdir, [output_prefix 'metadata-' num2str(ref_channel) '.mat']));
            
            stitch_time_slice(indir, tilefiles,...
                outdir, output_prefix,...
                ref_channel, NaN, NaN, 'Max', 0,... 
                save_stitched_image, assemble_from_metadata,...
                fullfile(outdir, ['log_' output_prefix num2str(c) '.txt']),...
                10, 10, c);
        end
        
    else
        for c = [ref_channel, setdiff(1:nChannels, ref_channel)]
            if c == ref_channel
                assemble_from_metadata = 0;
            else
                assemble_from_metadata = 1;
            end
            
            stitch_time_slice(indir, tilefiles,...
                outdir, output_prefix,...
                ref_channel, 7, 2, 'Max', 0,...
                save_stitched_image, assemble_from_metadata,...
                fullfile(outdir, ['log_' output_prefix num2str(c) '.txt']),...
                10, 10, c);
            
        end
    end
end
end

