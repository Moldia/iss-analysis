% iss base color reading of already aligned images
% alignment is done separately
% Xiaoyan, 2017


%% process reference anchor stain image
Iref = double(imread('b1_c2.tif'))/65535;

% tophat
Itop = imtophat(Iref, strel('disk', 3));

% binarize
Ibw = im2bw(Itop, .0025);

% declustering
Iws = watershed(-Iref);
Iws = double(Iws) .* double(Ibw);
Ibw = logical(Iws);

% remove any blob that is smaller than 3 px in area
Ibw = bwareaopen(Ibw, 3);
Label = bwlabel(Ibw);

% basic properties of reference blobs
refprops = regionprops(Label, 'Centroid', 'Area');
xypos = cat(1, refprops.Centroid);
nBlobs = length(refprops);

% visualization
figure, imshow(Iref, []);
hold on;
plot(xypos(:,1), xypos(:,2), 'r+');

%% read base intensities
filename_prefix = {};

intensities = zeros(nBlobs, 7, 4);

% read max intensity from individual objects
for b = 1:4     % base
    
    % first column is base number
    intensities(:,1,b) = b;
    % second column is object number
    intensities(:,2,b) = reshape(1:nBlobs, [], 1);

    for c = 2:6     % channel (skip DAPI)
        % load aligned images
        im = double(imread(['aligned_b' num2str(b) '_c' num2str(c) '.tif']))/65535;
        
        % tophat
        im = imtophat(im, strel('disk', 3));
        
        % use reference blob label image to measure intensities of current
        % base image
        baseprops = regionprops(Label, im, 'MaxIntensity');
        
        % store intensity data
        intensities(:,c+1,b) = cat(1,baseprops.MaxIntensity);
    end
end

intensities = permute(intensities, [1,3,2]);
intensities = reshape(intensities, [], 7);

%% write intensity reading file
fid = fopen('blob_intensities.csv', 'w');
fprintf(fid, lineformat('%s', 7), [{'base', 'object_num'}, catstrnum('channel', 2:6)]);
fprintf(fid, lineformat('%d', 7), intensities')
fclose(fid);



