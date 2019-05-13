%% Tile one image series
%  Xiaoyan, 2015
%  deprecated


%% parameters
imagename_prefix = 'hyb2_s1c';
imagename_suffix = '_ORG.tif';
series_num = 5;
Tilesize_x = 2000;
Tilesize_y = 2000;

%% read and pad the image
for s = 1:series_num
    I = imread([imagename_prefix, num2str(s), imagename_suffix]);
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
    
    % tile the image
    output_dir = ['Tiled_',strtok([imagename_prefix, num2str(s), imagename_suffix],'.')];
    mkdir(output_dir);
    for i = 1:tilenum(1)
        for j = 1:tilenum(2)
            Itemp = I((i-1)*Tilesize_y+1:i*Tilesize_y,(j-1)*Tilesize_x+1:j*Tilesize_x,:);
            imwrite(Itemp,[output_dir '\Tile' num2str((i-1)*tilenum(2)+j) '.tif'],'tif');
        end
    end
    
end