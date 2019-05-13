function [ROI_count,ROI_freq,ROI_proportion,ROI_area,In] = ...
    countinROI_f(ROI_number,Coord,Pos,name_uni,re_idx,p)
% count transcripts within ROI polygons
% Xiaoyan, 2015-8-26


ROI_count=[];
ROI_freq=[];
ROI_proportion=[];
ROI_area=[];
In=[];

Pos = correctcoordinates_f(Pos,1);

for i = 1:ROI_number
    poly = Coord(2:3,i);
    poly = ([poly{1};poly{2}])';
    in = inpolygon(Pos(:,1),Pos(:,2),poly(:,1),poly(:,2));
    
    In = [In,in];
    
    [m,n] = hist(re_idx(in),1:length(name_uni));
    ROI_count = [ROI_count,m'];
    ROI_freq = [ROI_freq,m'/sum(m)];
    ROI_proportion = [ROI_proportion,(m./p)'];
    ROI_area = [ROI_area,polyarea(poly(:,1),poly(:,2))];
end
In = logical(sum(In,2));    % used in plotting