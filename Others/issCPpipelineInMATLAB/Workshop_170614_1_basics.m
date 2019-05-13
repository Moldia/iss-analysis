% ISS image analysis workshop, 170614
% Xiaoyan
% tested on MATLAB R2016b


%% example images
% MIP image subset
% Imip = cell(4,6);
% for b = 1:4
%     imfile = ['I:\MIP\4028\', '161101_4028_36-1_b', num2str(b), '_mip.czi'];
%     
%     bfreader = loci.formats.Memoizer(bfGetReader(), 0);
%     % Initialize the reader with an input file
%     bfreader.setId(imfile);
%     
%     % get tile 89
%     bfreader.setSeries(88);
%     for c = 1:6
%         iPlane = bfreader.getIndex(0, c-1, 0)+1;
%         temp = bfGetPlane(bfreader, iPlane); 
%         Imip{b,c} = temp(1500:2000, 1500:2000);
%         imshow(Imip{b,c},[]);
% %         pause
%     end
%     
%     bfreader.close();
% end
% save WS170614 Imip


%% dummy image tophat
Ichips = imread('coloredChips.png');
Ichips = Ichips(1:200,191:390,:);
Ichips_gray = double(rgb2gray(Ichips))/255;
Ichips_bw = Ichips_gray<graythresh(Ichips_gray);
Ichips_bw = imfill(Ichips_bw, 'holes');

% randomly add 50 noise dots
Ichips_bw(randi(40000, 50, 1)) = 1;

% figure
figure(1); clf;
set(gcf, 'name', 'Binary dummy image', 'units', 'normalized', 'position', [.1 .2 .8 .6]); 

% binary image
subplot(231); imshow(Ichips_bw); title('original binary image');
axis on
pause()

% kernel
se = strel('disk', 9);
subplot(232); imshow(se.Neighborhood, []);
axis on
title('structuring element (kernel)');
pause()
set(gca, 'xlim', [-96 103], 'ylim', [-96 103], 'color', 'k');
pause()

% eroded
subplot(233);
imshow(imerode(Ichips_bw, se), []); title('erosion');
axis on
pause()

% dilated
subplot(234);
imshow(imdilate(Ichips_bw, se), []); title('dilation');
axis on
pause()

% opened
subplot(235);
Ichips_open = imopen(Ichips_bw, se);
imshow(Ichips_open, []); title({'opening'; 'erosion -> dilation'});
axis on
pause()

% closed
subplot(236);
imshow(imclose(Ichips_bw, se), []); title({'closing'; 'dilation -> erosion'});
axis on

%% base1 visulization
load WS170614 Imip
figure(1);
set(gcf, 'name', 'Base1 MIP images', 'units', 'normalized', 'position', [.1 .2 .8 .6]);
channels = {'DAPI', 'anchor', 'T', 'G', 'C' 'A'};
Ax = [];
for c = 1:6
    ax = subplot(2,3,c);
    Ax = [Ax, ax];
    imshow(Imip{1,c}, []);
    colorbar
    title(channels{c});
end
linkaxes(Ax, 'xy');

%% blob image tophat
figure(1); clf;
set(gcf, 'name', 'Base1 anchor', 'units', 'normalized', 'position', [.1 .2 .8 .6]); 

% original image
subplot(231); imshow(Imip{1,2}*30); title('anchor MIP image');
pause()
I = Imip{1,2}(301:330,301:330); colorbar
imshow(I, []); title('anchor MIP image');
axis on
colorbar
pause()

% kernel
se = strel('disk', 2);
subplot(232); imshow(se.Neighborhood, []);
axis on
title('structuring element (kernel)');
pause()
set(gca, 'xlim', [-11 18], 'ylim', [-11 18], 'color', 'k');
pause()

% eroded
subplot(233);
imshow(imerode(I, se), []); title('erosion'); colorbar;
axis on

% dilated
subplot(234);
imshow(imdilate(I, se), []); title('dilation'); colorbar;
axis on

% opened
subplot(235);
imshow(imopen(I, se), []); title({'opening'; 'erosion -> dilation'}); colorbar;
axis on

% tophatted
pause()
subplot(236);
imshow(imtophat(I, se), []); title({'tophat'; 'original - opening'}); colorbar;
axis on

%% structuring elements ~ tophat
figure(2);
set(gcf, 'name', 'SE and tophat', 'units', 'normalized', 'position', [.12 .22 .8 .6]); 

% Texas Red image (autofluorescence)
Ax = [];
ax = subplot(231); imshow(Imip{1,5}, []); title('TexRed MIP image');
Ax = [Ax, ax];
axis on
colorbar

% kernel
subplot(232);
s = se.Neighborhood;
s = padimg(s, 12, 12, 'NEWS');
bh = bar3(s, 1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
%     bh(i).EdgeColor = 'none';
end
colormap summer
title('"tophat"')

% tophatted
ax = subplot(233); 
Ax = [Ax, ax];
imshow(imtophat(Imip{1,5}, se), []); title({'tophat'; 'original - opening'}); colorbar;
axis on
linkaxes(Ax, 'xy');
pause()

% another shape - removing certain features
Ax2 = [];
Itombs = imread('2038248-690316_20130502_004.jpg');
Itombs = Itombs(:,400:end,:); 
ax = subplot(234); imshow(Itombs, []); title('random image');
Ax2 = [Ax2, ax];
axis on

SE = strel('rectangle', [15 10]);
subplot(235);
s = SE.Neighborhood;
s = padimg(s, 20, 20, 'NEWS');
bh = bar3(s, 1);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
    bh(i).EdgeColor = 'none';
end
colormap summer
title('"tophat"')

ax = subplot(236);
imshow(cat(3, imtophat(Itombs(:,:,1), SE), imtophat(Itombs(:,:,2), SE),imtophat(Itombs(:,:,3), SE)));
title({'tophat'; 'original - opening'});
Ax2 = [Ax2, ax];
axis on
linkaxes(Ax2, 'xy');


%% edges
figure(1); clf;
set(gcf, 'name', 'Edge detection', 'units', 'normalized', 'position', [.1 .2 .8 .6]); 

% binary image
subplot(231); imshow(Ichips_open); title('opened binary image');

% edges from binary
subplot(232); imshow(Ichips_open - imerode(Ichips_open, strel('disk',1)), []);
title({'edges' 'original - eroded'});
pause()

% grayscale image
subplot(234); imshow(Ichips_gray); title('original grayscale image');
pause()

% derivatives based methods
subplot(235);
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(Ichips_gray, hy, 'replicate');
Ix = imfilter(Ichips_gray, hx, 'replicate');
imshow(sqrt(Ix.^2 + Iy.^2));
title('sobel operator')
pause()

% % kernels
% figure; 
% subplot(121); imshow(hx, []); colorbar;
% subplot(122); imshow(hy, []); colorbar;
% colormap gray
% close gcf

imshow(edge(Ichips_gray, 'Sobel')); title('thresholded Sobel');
pause()
imshow(edge(Ichips_gray, 'Prewitt')); title('thresholded Prewitt');
pause()
imshow(edge(Ichips_gray, 'Roberts')); title('thresholded Roberts');
pause()

% Laplacian of Gaussian
fh = fspecial('LoG', 25, 2);
subplot(236); imshow(imfilter(Ichips_gray, fh),[]);
colorbar
title('LoG filtered')
pause()

% kernel
subplot(233); surf(fh); title('LoG kernel');  pause()

subplot(236)
imshow(edge(Ichips_gray, 'log')); title('thresholded LoG');

%% watershed
figure(1); clf;
set(gcf, 'name', 'Watershed on intensity', 'units', 'normalized', 'position', [.1 .2 .8 .6]); 

Ax = [];

% base1 DAPI image
ax = subplot(231); Ax = [Ax, ax]; linkaxes(Ax, 'xy');
I = double(Imip{1,1}(1:200, 301:500))/65535;
imshow(I, []); title('base1 DAPI image');
pause()

% Otsu threshold
ax = subplot(232);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Ibw = im2bw(I, graythresh(I));
imshow(Ibw, []); title('Otsu thresholded');
pause()

% invert
ax = subplot(233);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
imshow(-I, []); title('inverted image');
pause()

% watershed
ax = subplot(234);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Iws = watershed(-I);
imshow(Iws, []); title('watershed on inverted');
pause()

% impose foreground
ax = subplot(235);   Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Iws = double(Iws) .* double(Ibw);
imshow(label2rgb(Iws), []); title('thresholded imposed on watershed');
pause()

% h-minima before watershed
ax = subplot(236); Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Ismooth = imfilter(I, fspecial('gaussian', 20));
Ihmin = imhmin(-Ismooth, .01);
imshow(Ihmin, []); title('inverted, h-min transformed');
pause()

% watershed
subplot(234);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Iws2 = watershed(Ihmin);
imshow(Iws2, []); title('watershed on inverted 2');
pause()

% impose binary
subplot(235);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Iws2 = double(Iws2) .* double(Ibw);
imshow(label2rgb(Iws2), []); title('thresholded imposed on watershed 2');

%% distance transform and watershed
figure(1); clf;
set(gcf, 'name', 'Watershed on distance transformed', 'units', 'normalized', 'position', [.1 .2 .8 .6]); 

Ax = [];

% grayscale image
ax = subplot(231); Ax = [Ax, ax]; linkaxes(Ax, 'xy');
I = double(Imip{1,1}(1:200, 301:500))/65535;
I = imtophat(I, strel('disk', 20));
imshow(I, []); title('grayscale image');

% thresholding
ax = subplot(232);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Ibw = im2bw(I, graythresh(I));
imshow(Ibw, []); title('Otsu thresholded');
pause()

% distance transform
ax = subplot(233);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Idist = bwdist(~Ibw);
imshow(Idist, []); title('distance transform'); drawnow
pause(3)
imshow(-Idist, []); title('inverted distance transform'); drawnow
pause()

% watershed
ax = subplot(234);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Iws = watershed(-Idist);
imshow(Iws, []); title('watershed on distance'); drawnow
pause(3)

% impose binary
Iws = double(Iws) .* double(Ibw);
imshow(label2rgb(Iws), []); title('thresholded imposed on watershed'); drawnow
pause()

% force seed positions based on local intensity maxima
ax = subplot(235);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
maximaDapi = ones(size(Idist));
maxfilter = imdilate(Ismooth, strel('disk', 9));
maximaDapi(Ismooth < maxfilter) = 0;
maximaDapi = imdilate(maximaDapi, ones(2));
imshow(maximaDapi, []); title('forced seeds'); drawnow
pause()

imposedDapi = imimposemin(-Idist, maximaDapi);   % impose minima, force seed positions
imshow(imposedDapi, []); title('"landscape"'); drawnow
pause()

ax = subplot(236);  Ax = [Ax, ax]; linkaxes(Ax, 'xy');
Iws2 = watershed(imposedDapi);
Iws2 = double(Iws2) .* double(Ibw);
imshow(label2rgb(Iws2), []); title('watershed after imposed minima'); drawnow

% force seed positions based on distance transform
subplot(235);
maximaDapi = ones(size(Idist));
maxfilter = imdilate(Idist, strel('disk', 5));
maximaDapi(Idist < maxfilter) = 0;
maximaDapi = imdilate(maximaDapi, ones(2));
imshow(maximaDapi, []); title('forced seeds'); drawnow
pause()

imposedDapi = imimposemin(-Idist, maximaDapi);   % impose minima, force seed positions
imshow(imposedDapi, []); title('"landscape"'); drawnow
pause()

subplot(236);
Iws2 = watershed(imposedDapi);
Iws2 = double(Iws2) .* double(Ibw);
imshow(label2rgb(Iws2), []); title('watershed after imposed'); drawnow
