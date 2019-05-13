%% stitch images
% default output: 8bit grayscale or 24bit RGB
% Xiaoyan, 2016-2-19
% deprecated, new version as function

%% input
tilesize = 2000;
ntilesX = 15;
ntilesY = 11;
newScale = 0.5;

name_prefix = 'E:\PROOOJECTS\12_Neuron_mapping\150616_FF\4028_1_2\CP_150702_Cell\Outlines\CellOutlines';   % include path
name_suffix = '.png';   % include file format
digit_num = 3;  % total digits in the file name, eg. 4 if image0001.tif; enter 0 if eg. image1.tif, ..., image 10.tif
output_name = 'test.png';   % include path, name and format    

%% start
t = 1;
tile = '1';
while t < digit_num
    tile = ['0', tile];
    t = t+1;
end
imtemp = imread(strcat(name_prefix,tile,name_suffix));
if strcmp(class(imtemp),'uint8') 
    bit8 = 1;
else
    bit8 = 0;
end
d = size(imtemp);
if length(d) == 3
    finalimage = zeros(ntilesY*tilesize*newScale,ntilesX*tilesize*newScale,3,'uint8');
else
    finalimage = zeros(ntilesY*tilesize*newScale,ntilesX*tilesize*newScale,'uint8');
end

for i = 1:ntilesY
    for j = 1:ntilesX
        t = (i-1)*ntilesX+j;
        if t ~= 1
            tile = num2str(t);
            t = floor(log10(t))+1;
            while t < digit_num
                tile = ['0', tile];
                t = t+1;
            end
            imtemp = imread(strcat(name_prefix,tile,name_suffix));
        end
        if ~bit8
            imtemp = uint8(imtemp/65535*255);
        end
        imtemp = imresize(imtemp,newScale);
        
        finalimage((i-1)*tilesize*newScale +1:...
            i*tilesize*newScale,...
            (j-1)*tilesize*newScale+1:...
            j*tilesize*newScale,:) = imtemp;
    end
end

imshow(finalimage);
imwrite(finalimage,output_name);

