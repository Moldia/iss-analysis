function [CellMap, CellYX] = iss_segment_dapi(o, Dapi0)
% CellMap = iss_segment_dapi(DapiIm)
%
% segments a DAPI image and assigns each pixel to a cell. 
%
% Output CellMap is same size as input DapiIm, with integer entries for
% each pixel, saying which cell it belongs to. (Zero if unassigned)
%
% Output CellsYX is YX coordinates of each cell's center (Ncells by 2)


%Dapi = imadjust(imfilter(Dapi0, fspecial('gaussian', 3))); % filter and contrast enhancement
Dapi = imadjust(Dapi0);
ImSz = size(Dapi);
Debug = 0;
%% threshold the map
ThreshVal = prctile(Dapi(:), o.DapiThresh);

%bwDapi = imerode(imfill(Dapi>ThreshVal, 'holes'), strel('disk', 2));
bwDapi = imerode(Dapi>ThreshVal, strel('disk', 2));

if Debug
    figure(300)
    subplot(2,1,1);
    imagesc(Dapi); 
    subplot(2,1,2);
    imagesc(bwDapi);
    colormap bone
end
%% find local maxima 
dist = bwdist(~bwDapi);
dist0 = dist;
dist0(dist<o.DapiMinSize)=0;
ddist = imdilate(dist0, strel('disk', o.DapiMinSep));
% clear dist dist0
dLocMaxIm = imregionalmax(ddist); % with dilation


impim = imimposemin(-dist0, dLocMaxIm);

if Debug
    figure(301);
    subplot(2,2,1)
    imagesc(dist);

    subplot(2,2,2)
    subplot(2,2,3)
    imagesc(dLocMaxIm);
    
    subplot(2,2,4)
    imagesc(impim);
end
%% segment
% remove pixels at watershed boundaries
bwDapi0 = bwDapi;
bwDapi0(watershed(impim)==0)=0;

% assign all pixels a label
labels = uint32(bwlabel(bwDapi0));
[d, idx] = bwdist(bwDapi0);

% get rid of too small regions
% rprops = regionprops(labels);
% DeleteMe = [rprops.Area]<o.MinDapiArea;


% now expand the regions by a margin
CellMap0 = zeros(ImSz, 'uint32');
Expansions = (d<o.DapiMargin);
CellMap0(Expansions) = labels(idx(Expansions));

rProps = regionprops(CellMap0); % returns XY coordinate and area
BigEnough = [rProps.Area]>o.MinCellArea;
NewNumber = zeros(length(rProps),1);
NewNumber(~BigEnough) = 0;
NewNumber(BigEnough) = 1:sum(BigEnough);
CellMap = CellMap0;
CellMap(CellMap0>0) = NewNumber(CellMap0(CellMap0>0));

CellYX = fliplr(vertcat(rProps(BigEnough).Centroid)); % because XY


if Debug
     figure(302); clf
%     subplot(2,2,1);
%     image(label2rgb(labels, 'jet', 'w', 'shuffle'));
% 
%     subplot(2,2,2);
    imagesc(Dapi); colormap(bone); hold on
    plot(CellYX(:,2), CellYX(:,1), 'rx');

%     subplot(2,2,3);
%     imagesc(labels);
% 
%     % now give every pixel to its nearest neighbor
%     subplot(2,2,4);
%     image(CellMap);
end
%% make image with dashed boundaries
Boundaries = (CellMap ~= imdilate(CellMap,strel('disk', 1)));
DapiBoundaries = Dapi;

OddPix = mod((1:size(Dapi,2)) + (1:size(Dapi,1))', 2);
DapiBoundaries(Boundaries & OddPix) = .3 * max(Dapi(:));
DapiBoundaries(Boundaries & ~OddPix) = 0;


imwrite(DapiBoundaries, fullfile(o.OutputDirectory, 'background_boundaries.tif'));
save CellMap CellMap
% figure(904375)
% imagesc(DapiBoundaries)


return

%% make a pseudocolored image

colors = label2rgb(CellMap, 'hsv', 'w', 'shuffle');

NormDapi = min(double(Dapi)/double(max(Dapi(:))),1);
pColorIm = bsxfun(@times, double(colors)/255, NormDapi);
Boundaries = (CellMap ~= imdilate(CellMap,strel('disk', 1)));
[y, x] = find(Boundaries);
for i=1:3
    pColorIm(sub2ind(size(pColorIm), y, x, i*ones(size(y))))=.3; %ones(size(y))*i]))=1;
end

figure(904375)
image(pColorIm)

%imdilate(bwDapiSep, strel('disk', 1))-bwDapiSep;
%pColorIm

return
% subplot(2,2,3)
% Props = regionprops(bwDapiSep, 'EquivDiameter');
% histogram([Props.EquivDiameter],0:max([Props.EquivDiameter]));
% 
% 
% excludeDapi = cellfun(@(v) sum(maximaDapi(v))==0, propDapi(2,:));       % no maxima
% excludeDapi = excludeDapi | ~cellfun(@(v) v>8 & v<100, propDapi(1,:));  % size threshold
%%
%Thresh = 500;%prctile(Dapi(:), 90);
Thresh = graythresh(Dapi)*max(Dapi(:))
bim = (Dapi>Thresh);
imagesc(bim);
colormap bone
%%
figure(302)
bwdim = bwdist(~bim);
imagesc(bwdim);
colorbar
%%
% figure(303)
% Outside = find(Dapi<Thresh);
% gdim = graydist(Dapi, Outside);
% imagesc(gdim)
%%
figure(304)
imax = imregionalmax(bwdim);
imagesc(imax.*bwdim)
%%
figure(305)
bwdim2 = bwdim; 
bwdim2(imax) = bwdim2(imax)-1;
imax2 = imregionalmax(bwdim2);
imagesc(imax2.*bwdim2);
%% 
figure(306);
iimp = imimposemin(-bwdim, imax2 & (bwdim>9));
imagesc(iimp);
%%
L = watershed(iimp);
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
image(bsxfun(@times, double(Lrgb), double(bim)));

%%
se = strel('disk', 4);
%grad = imdilate(Dapi, se) - imerode(Dapi, se);
hy = fspecial('sobel');
gy = double(imfilter(Dapi, hy, 'replicate'));
gx = double(imfilter(Dapi, hy', 'replicate'));
grad = sqrt(gy.^2 + gx.^2);
imagesc(grad)