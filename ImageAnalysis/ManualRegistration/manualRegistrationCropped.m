%TODO: read nuclei image
base=imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\DO\25_c1.tif');
%Crop the base image and save
base=base(70:size(base,1), 20:size(base,2)-119);
imwrite(base, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\Tranformed\DO\25_c1.tif');

%Base image
base=imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\DO\25_c2.tif');
%Crop the base image and save
base=base(70:size(base,1), 20:size(base,2)-119);
imwrite(base, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\Tranformed\DO\25_c2.tif');

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\1\25.1_c3.tif');
input=input(70:size(base,1), 20:size(base,2)-119); %Crop
[pointsInput, pointsBase] = cpselect(input, base, 'Wait', true);
t = cp2tform(pointsInput, pointsBase, 'affine');
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\Tranformed\1\25.1_c3.tif');
rgbim = cat(3, base, transformed, zeros(size(base)));
figure, imshow(rgbim)

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\1\25.1_c1.tif');
input=input(70:size(base,1), 20:size(base,2)-119);
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\Tranformed\1\25.1_c1.tif');

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\1\25.1_c2.tif');
input=input(70:size(base,1), 20:size(base,2)-119);
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\Tranformed\1\25.1_c2.tif');


input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\1\25.1_c4.tif');
input=input(70:size(base,1), 20:size(base,2)-119);
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\25\Tranformed\1\25.1_c4.tif');