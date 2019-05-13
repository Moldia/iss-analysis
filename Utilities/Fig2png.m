% save fig file as raster image
% Xiaoyan, 2017

markersize = 10;

% get child objects of current axis
ch = get(gca, 'children');

% number of line objects
iMarker = [];
for i = 1:length(ch)
    if strcmp(ch(i).Type, 'line') && ~isempty(ch(i).DisplayName) && strcmp(ch(i).Visible, 'on')
        iMarker = [iMarker, i];
    end
end

% get image and convert to color 
im = ch(end);
im = im.CData;
if numel(size(im)) == 2     % 2D image
    im = repmat(im,1,1,3);
end

% get original color
col = [];
for i = 1:length(iMarker)
    col(i,:) = ch(iMarker(i)).Color;
end

% % process background image
% blobimage = im(:,:,1);
% blobimage = imtophat(blobimage, strel('disk', 7));
% im(:,:,1) = blobimage*8;
% figure, imshow(im);

% save original spots as a single pixel object
for i = 1:length(iMarker)
    temp = ch(iMarker(i));
    for j = 1:length(temp.XData)
        for k = 1:3
            im(round(temp.YData(j))-ceil(markersize/2):round(temp.YData(j))+ceil(markersize/2),...
                round(temp.XData(j))-ceil(markersize/2):round(temp.XData(j))+ceil(markersize/2),k) = (255^(1+isa(im, 'uint16')))*col(i,k);
        end
    end
end

figure, imshow(im);
imwrite(im, 'OriginalResolution.tif');
