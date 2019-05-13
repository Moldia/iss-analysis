%% Pad and tile a single image
%  Xiaoyan, 2015-6-10
%  deprecated


%% parameters
image = 'base2_c6_ORG.tif';
Tilesize_x_original = 2000;
Tilesize_y_original = 2000;
output_directory = 'Tiled_c6\';

%% read the image
I = imread(image);
Isize = size(I);
if length(Isize)==3
    col = 1;
    Isize = Isize(1:2);
else
    col = 0;
end

%% pad the left upper part of the image
if col
    I = [zeros(Tilesize_y_original/2,Isize(2),3);I];
    I = [zeros(Isize(1)+Tilesize_y_original/2,Tilesize_x_original/2,3),I];  
    Isize = size(I);
    Isize = Isize(1:2);
else
    I = [zeros(Tilesize_y_original/2,Isize(2));I];
    I = [zeros(Isize(1)+Tilesize_y_original/2,Tilesize_x_original/2),I];
    Isize = size(I);
end

%% pad and tile as in normal tiling
tilenum = ceil(Isize./[Tilesize_y_original,Tilesize_x_original]);
padsize = tilenum.*[Tilesize_y_original,Tilesize_x_original]-Isize;
if col
    I = [I;zeros(padsize(1),Isize(2),3)];
    I = [I,zeros(size(I,1),padsize(2),3)];
else
    I = [I;zeros(padsize(1),Isize(2))];
    I = [I,zeros(size(I,1),padsize(2))];
end

%% tile the image
if strcmp(output_directory(end),'\')
    output_directory = output_directory(1:end-1);
end
mkdir(output_directory);
for i = 1:tilenum(1)
    for j = 1:tilenum(2)
        Itemp = I((i-1)*Tilesize_y_original+1:i*Tilesize_y_original,(j-1)*Tilesize_x_original+1:j*Tilesize_x_original,:);
        imwrite(Itemp,[output_directory '\Tile' num2str((i-1)*tilenum(2)+j) '.tif'],'tif');
    end
end