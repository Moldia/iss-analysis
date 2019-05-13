% create image subset containing ROI and extract reads within polygon
% Xiaoyan, 2017

polygon = csvread('polygon_coordinates.csv', 1);
Idapi = imread('Ab_c1_Aligned.tif');
Iab = imread('Ab_c3_Aligned.tif');

[name, pos] = getinsitudata('SpotCoordinates.csv', 1, 1);
roundedCoord = round(pos);

for i = 1:max(polygon(:,1))
    coordinates = polygon(polygon(:,1)==i,2:3);
    
    % a bit of buffer, can be excluded, used to create an image subset
    % containing only ROI+-100 px
    imin = floor(min(coordinates(:,2))-100);
    imax = ceil(max(coordinates(:,2))+100);
    jmin = floor(min(coordinates(:,1))-100);
    jmax = ceil(max(coordinates(:,1))+100);
    
    mask = poly2mask(coordinates(:,1)-jmin, coordinates(:,2)-imin, imax-imin+1, jmax-jmin+1);
    
    % square area
    insquare = readsinsqr(pos, [jmin imin jmax imax]);
    subpos = bsxfun(@minus, pos(insquare,:), [jmin, imin]);
    insquare = find(insquare);

    % polygon
    subpos(subpos<.5) = .5;
    inroi = logical(readsinroi(subpos, mask));
    inroi = insquare(inroi);
    
    imshow(mask); hold on; plot(pos(inroi,1)-jmin, pos(inroi,2)-imin, '.');
    
    % save subset ROI images
    Dapi = Idapi(imin:imax,jmin:jmax).*uint16(mask);
    Ab = Iab(imin:imax,jmin:jmax).*uint16(mask);
    imwrite(Dapi, ['Ab_c1_ROI' num2str(i) '.tif']);
    imwrite(Ab, ['Ab_c3_ROI' num2str(i) '.tif']);
    
    % write gene name and coordinate of reads in polygon, and original
    % index (coordinate system matching subset ROI image)
    towrite = [name(inroi), num2cell([pos(inroi,1)-jmin, pos(inroi,2)-imin, inroi])]';
    fid = fopen(['spots_ROI' num2str(i) '.csv'], 'w');
    fprintf(fid, 'name, sub_x, sub_y, original_index\n');
    fprintf(fid, '%s,%d,%d,%d\n', towrite{:});
    fclose(fid);
    
end