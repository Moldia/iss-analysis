function [ROI_number,Coord,blank] = getcsvROIinfo_f(coordinate_file)
% get ROI coordinates from a csv file
% example file:
% Polygon id,x coordiates,y coordinates
% 2,1.153496e+04,1.904365e+03
% 2,1.009527e+04,2.444246e+03
% 2,8.595605e+03,4.993684e+03
% 2,1.021525e+04,5.233631e+03
% 2,1.213482e+04,2.654200e+03
% 3,1.318459e+04,1.633118e+04
% 3,1.021525e+04,1.642116e+04
% 3,1.177490e+04,1.873065e+04
% 3,1.363449e+04,1.828075e+04
% 4,1.036521e+04,1.369177e+04
% 
% Xiaoyan, 2014-12-8

Coord_file = csvread(coordinate_file,1);
ROI_number = max(Coord_file(:,1));
Coord = {}; blank = [];
j = 0;
for i = 1:ROI_number
    if ~isempty(find(Coord_file(:,1)==i))
        if length(find(Coord_file(:,1)==i))==1 || length(find(Coord_file(:,1)==i))==2
            blank = [blank;i];
        else
            j = j+1;
            Coord{1,j} = ['roi' num2str(i)];
            Coord{2,j} = (Coord_file(Coord_file(:,1)==i,2))';
            Coord{2,j} = [Coord{2,j},Coord{2,j}(1)];    % connect the first and last point
            Coord{3,j} = (Coord_file(Coord_file(:,1)==i,3))';
            Coord{3,j} = [Coord{3,j},Coord{3,j}(1)];    % connect the first and last point
        end
    else
        blank = [blank;i];
    end
end
