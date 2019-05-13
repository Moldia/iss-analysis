function [ROI_number,Coord] = getIJROIinfo_f(coordinate_file)
% get ROI coordinates from a text file copied from Image J
% example file:
% roi1
% 15864,18384,19272,19032,17808,16128,14208,13032,12720,13752,15720
% 18360,18576,20208,22440,23472,23688,23568,21888,20160,19104,18360
% 
% necrotic
% 8784,8712,9408,10560,11328,11376,9504
% 17640,16152,16128,16752,17472,18048,18096
% 
% Xiaoyan, 2014-11-27

fid = fopen(coordinate_file);
Tcoordinate = [];
textline = fgets(fid);
while textline ~= -1
    if isempty(textline)
    else
        Tcoordinate = [Tcoordinate;{textline}];
    end
    textline = fgets(fid);
end
fclose(fid);

ROI_number = (length(Tcoordinate)+1)/4;
Coord = cell(3,ROI_number);
for i = 1:ROI_number
    Coord(1,i) = Tcoordinate(4*(i-1)+1);
    Coord(1,i) = strtok(Coord(1,i));
    Coord(2,i) = {cellfun(@str2num,strsplit(Tcoordinate{4*(i-1)+2},','))};
    Coord{2,i} = [Coord{2,i},Coord{2,i}(1)];    % connect the first and last point
    Coord(3,i) = {cellfun(@str2num,strsplit(Tcoordinate{4*(i-1)+3},','))};
    Coord{3,i} = [Coord{3,i},Coord{3,i}(1)];    % connect the first and last point
end