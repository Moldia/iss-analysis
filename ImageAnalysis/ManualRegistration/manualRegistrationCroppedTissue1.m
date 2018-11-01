%Base image
base=imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\Tranformed\DO\1_c2.tif');

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\2\1.2_c4.tif');
input=input(70:size(base,1), 20:size(base,2)-119); %Crop
[pointsInput, pointsBase] = cpselect(input, base, 'Wait', true);
t = cp2tform(pointsInput, pointsBase, 'affine');
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\Tranformed\2\1.2_c4.tif');
rgbim = cat(3, base, transformed, zeros(size(base)));
figure, imshow(rgbim)

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\2\1.2_c1.tif');
input=input(70:size(base,1), 20:size(base,2)-119);
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\Tranformed\2\1.2_c1.tif');

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\2\1.2_c2.tif');
input=input(70:size(base,1), 20:size(base,2)-119);
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\Tranformed\2\1.2_c2.tif');


input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\2\1.2_c3.tif');
input=input(70:size(base,1), 20:size(base,2)-119);
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Tissue manual alignment\12392\1\Tranformed\2\1.2_c3.tif');