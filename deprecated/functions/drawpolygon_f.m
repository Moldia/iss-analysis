function [Coord,Coord_write,blank] = drawpolygon_f(I,ROI_number,varargin)
% Draw polygons from an image and get coordinates
% Xiaoyan, 2014-12-15

Coord = {}; Coord_write = []; blank = [];

switch length(varargin)
    case 0
        scale = 1;
    case 1
        scale = varargin{1};
end
    
for i = 1:ROI_number
    figure;
    set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.85]);
    
    imshow(I,[]);
    
%         h = imfreehand;
    h = impoly;
    pos = getPosition(h);
    % wait until double click to confirm polygon
    position = wait(h);
    position = floor(position/scale);
    Coord = [Coord,...
        [{['ROI' num2str(i)]};...
        {[position(:,1);position(1,1)]'};...
        {[position(:,2);position(1,2)]'}]];
    Coord_write = [Coord_write; zeros(size(position,1),1)+i,position];
    if size(position,1)==1 || size(position,1)==2
        blank = [blank;i];
    end
    close(gcf);
end