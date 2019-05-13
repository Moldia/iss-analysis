%% Tile a single image
%  Xiaoyan, 2015-4-18
%  deprecated

%% parameters
image = 'Kseq2_HE.TIF';
Tilesize_x = 2000;
Tilesize_y = 3000;
output_directory = 'Tiled\';

%% read and pad the image
I = imread(image);
Isize = size(I);
if length(Isize)==3
    col = 1;
    Isize = Isize(1:2);
else
    col = 0;
end
tilenum = ceil(Isize./[Tilesize_y,Tilesize_x]);
padsize = tilenum.*[Tilesize_y,Tilesize_x]-Isize;
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
        Itemp = I((i-1)*Tilesize_y+1:i*Tilesize_y,(j-1)*Tilesize_x+1:j*Tilesize_x,:);
        imwrite(Itemp,[output_directory '\Tile' num2str((i-1)*tilenum(2)+j) '.tif'],'tif');
    end
end