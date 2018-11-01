base=imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\DO\1.DO_c3m05_ORG.tif');
input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c2m05_ORG.tif');
[pointsInput, pointsBase] = cpselect(input, base, 'Wait', true);
t = cp2tform(pointsInput, pointsBase, 'affine');
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c2m05_ORG.tif');
rgbim = cat(3, base, transformed, zeros(size(base)));
figure, imshow(rgbim)

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c3m05_ORG.tif');
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c3m05_ORG.tif');

input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c4m05_ORG.tif');
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c4m05_ORG.tif');


input = imread('C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c1m05_ORG.tif');
transformed = imtransform(input, t, 'XData', [1 size(base, 2)], 'YData', [1 size(base,1)]);
imwrite(transformed, 'C:\Users\marmi748.RUDBECKLAB\Desktop\Cell lines manual alignment\LNCaP\1\5\1\1.1_c1m05_ORG.tif');