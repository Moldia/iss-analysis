% save fig file as raster image
% Xiaoyan, 2017


% get child objects of current axis
ch = get(gca, 'children');

% number of line objects
nMarkers = 0;
for i = 1:length(ch)
    if strcmp(ch(i).Type, 'line')
        nMarkers = nMarkers + 1;
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
for i = 1:nMarkers
    col(i,:) = ch(i).Color;
end

% % process background image
% blobimage = im(:,:,1);
% blobimage = imtophat(blobimage, strel('disk', 7));
% im(:,:,1) = blobimage*8;
% figure, imshow(im);

% save original spots as a single pixel object
for i = 1:nMarkers
    temp = ch(i);
    for j = 1:length(temp.XData)
        for k = 1:3
            im(round(temp.YData(j))-2:round(temp.YData(j))+2,...
                round(temp.XData(j))-2:round(temp.XData(j))+2,k) = 65535*col(i,k);
        end
    end
end

figure, imshow(im);
imwrite(im, '200_13_1.png');
