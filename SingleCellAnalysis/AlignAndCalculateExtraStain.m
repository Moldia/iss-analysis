% align Ab images to iss reference DAPI
% calculate mean Ab intensity of each cell object
% the cell objects are the same as used in single cell analysis
% Xiaoyan, 2018


%% tile with 10% buffer (2000 + 200 px)
I = imread('..\tiffs\base1_c1_ORG.tif');
refsize = size(I);
buffertile(I, 2000, 200, 'tiffs\BufferTiled\base1_c1_ORG');

for c = 1:3
    I = imread(['..\tiffs\AlignedImages\IF_c' num2str(c) '_ORG.tif']);
%     I2 = padimg(I, refsize(2)-size(I,2), -400, 'EN');
%     I2 = padimg(I2, 0, refsize(1)-size(I,1)+400, 'S');
    
%     if c == 1
%         clf; imshow(I2, []);
%         pause()    
%     end
%     
%     imwrite(I2, ['tiffs\IF-restain_c' num2str(c) '_ORG_refsize.tif']);

    % tile with 10% buffer each side
    buffertile(I, 2000, 200, ['..\tiffs\BufferTiled\IF-restain_c' num2str(c) '_ORG']);
end


%% tile alignment 
fls = cellstr(ls('..\tiffs\BufferTiled\base1_c1_ORG'));
fls = fls(cellfun(@(v) strncmp(v, 'tile', 4), fls));

mkdir('IF_aligned')

parfor i = 1:numel(fls) % with parallelization
%     clear ref tform
    
    ref = imread(['..\tiffs\BufferTiled\base1_c1_ORG\' fls{i}]);
    float = imread(['..\tiffs\BufferTiled\IF-restain_c1_ORG\' fls{i}]);
    
    % align with pixel resolution (translation only, should be fine for most cases)
    % fft and phase correlation
    tform = imregcorr(ref, float, 'translation');
    
    float_aligned = padimg(float, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
    
    % pad to original size (just for convenience)
    float_aligned = padimg(float_aligned,...
        size(ref,2)-size(float_aligned,2),...
        size(ref,1)-size(float_aligned,1), 'SE');
 
%     figure, imshowpair(ref, float_aligned);
    imwrite(cat(3, uint8(ref/255), uint8(float_aligned/255), zeros(size(ref), 'uint8'))*5,...
        ['IF_aligned\' strtok(fls{i}, '.') '.jpg']);

%     clear float float_aligned
    
    for c = 2:3
        float = imread(['..\tiffs\BufferTiled\IF-restain_c' num2str(c) '_ORG\' fls{i}]);
        
        % use the same transformation matrix
        float_aligned = padimg(float, -round(tform.T(3,1)), -round(tform.T(3,2)), 'NW');
        
        % pad
        float_aligned = padimg(float_aligned,...
            size(ref,2)-size(float_aligned,2),...
            size(ref,1)-size(float_aligned,1), 'SE');
        
        % write
        mkdir(['IF_aligned\IF-restain_c' num2str(c) '_ORG']);
        imwrite(float_aligned, ['IF_aligned\IF-restain_c' num2str(c) '_ORG\' fls{i}]);
        
%         clear float float_aligned
    end
    
end

%% stitch back red and green IF images
% ref = imread('tiffs\base1_c1_ORG.tif'); 
for c = 2:3
    im = stitch_buffered_tiles(...
        '..\tiffs\base1_c1_ORG.tif',...
        ['IF_aligned\IF-restain_c' num2str(c) '_ORG\tile'],...
        2000, 200, '.tif', 0);
    
%     figure, imshow(cat(3, ref*5, im*10, im*10))
    imwrite(im, ['IF_aligned\aligned_c' num2str(c) '.tif']);
end
    
%% load segmented cell label image from single cell analysis   
load ..\singleCellanalysis\Stitched\CellLookupTable.mat CellCorrectionTable
ntilesX = size(CellCorrectionTable, 2);
ntilesY = size(CellCorrectionTable, 1);

Ilabeled = cell(ntilesY,ntilesX);

% load tile nuclei label images
for i = 1:ntilesY
    fprintf('%.2f percet read.\n', (i-1)/ntilesY*100);
    parfor j = 1:ntilesX
        t = ntilesX*(i-1)+j;
        I1 = imread(['..\singleCellanalysis\CP_170825_Cell\Segmentation\NucleiLabel_', paddigits(t,3), '.tif']);
        I1 = relabelimg(I1);
        Ilabeled(i,j) = {I1};
    end
    
end  
clear I
fprintf('100 percet finished.\n');
% save('FullLabel.mat', 'Ilabeled', '-v7.3');

% relabel nuclei objects using correction table from single cell analysis
IRelabeled = cell(ntilesY,ntilesX);
for i = 1:ntilesY
    fprintf('%.2f percet finished.\n', (i-1)/ntilesY*100);
    parfor j = 1:ntilesX
        t = ntilesX*(i-1)+j;
        
        temp = uint32(Ilabeled{i,j});
        mapindex = CellCorrectionTable{i,j};
        
        for k = size(mapindex,1):-1:1
            temp(temp==mapindex(k,1)) = mapindex(k,2);
        end
        IRelabeled{i,j} = reshape(temp, 2000, []);

    end
    
end  

IRelabeled(cellfun(@isempty, IRelabeled)) = {zeros(2000, 'uint32')};
save('FullLabel.mat', 'IRelabeled', '-7.3');

%% stitch to full image and expand 20 px from nuclei
clear
load('FullLabel.mat', 'IRelabeled');
ntilesX = size(IRelabeled, 2);
ntilesY = size(IRelabeled, 1);

% stitch
Ifull = resizestitch([ntilesX, ntilesY], 2000, 1, IRelabeled, 1);
save('FullLabel.mat', 'Ifull', '-append');
clear IRelabeled

% distance transform
[D, idx] = bwdist(Ifull);
save('FullLabel.mat', 'D', 'idx', '-append');

% discard everything outside 20px expansion
D = D<=20;
L = zeros(size(Ifull), 'uint32');
L(D~=0) = Ifull(idx(D~=0));
save('FullLabel.mat', 'L', '-append');

% outlines
LOutlines = L~=imerode(L,strel('disk', 1));
save('FullLabel.mat', 'LOutlines', '-append');

clear Ifull idx D 

%% region properties
% first Ab image
I1 = imread('IF_aligned\aligned_c2.tif');

% cut edges in label image
L = padimg(L, size(I1,2)-size(L,2), size(I1,1)-size(L,1));
save('FullLabel.mat', 'L', '-append');

% calculate mean intensity of each cell
stats = regionprops(L, I1, 'MeanIntensity', 'MaxIntensity', 'Centroid');
save('FullLabel.mat', 'stats', '-append');

% second Ab image
I2 = imread('IF_aligned\aligned_c3.tif');
stats_2 = regionprops(L, I2, 'MeanIntensity', 'MaxIntensity');
save('FullLabel.mat', 'stats_2', '-append');

%% convert into simple arrays
centroid = cat(1, stats.Centroid);
meanIntensity_1 = cat(1, stats.MeanIntensity);
% maxIntensity = cat(1, stats.MaxIntensity);
meanIntensity_2 = cat(1, stats_2.MeanIntensity);

clear L stats stats_2

%% visualization
% first Ab
clf; 
imshow(imresize(I1, .5)*20, []);
hold on;
addimlayer(imresize(uint16(LOutlines)*65535, .5), .2);

% show top 5% mean intensity cells
toshow = meanIntensity_1>prctile(meanIntensity_1, 98);
text(centroid(toshow,1)*.5, centroid(toshow,2)*.5,...
    cellfun(@num2str, num2cell(round(meanIntensity_1(toshow), 0)), 'uni', 0),...
    'color', 'g', 'fontsize', 5);
title('channel 2, mean intensity top 2%')

% second Ab
clf; 
imshow(imresize(I2, .5)*20, []);
hold on;
addimlayer(imresize(uint16(LOutlines)*65535, .5), .2);

% show top 5% mean intensity cells
toshow = meanIntensity_2>prctile(meanIntensity_2, 95);
text(centroid(toshow,1)*.5, centroid(toshow,2)*.5,...
    cellfun(@num2str, num2cell(round(meanIntensity_2(toshow), 0)), 'uni', 0),...
    'color', 'r', 'fontsize', 5);
title('channel 3, mean intensity top 5%')

clear I1 I2

%% save to csv
data = importdata('..\singleCellanalysis\CellBlobs\CellBlobs_QT_0.35.csv');
towrite = [data.data, meanIntensity_1(data.data(:,1)), meanIntensity_2(data.data(:,1))];
columns = [data.colheaders, {'meanAbIntensityC2'}, {'meanAbIntensityC3'}];

fid = fopen('CellBlobs_QT_0.35_Ab.csv', 'w');
fprintf(fid, lineformat('%s', numel(columns)), columns{:});
fprintf(fid, lineformat('%d', numel(columns)), towrite');
fclose(fid);

