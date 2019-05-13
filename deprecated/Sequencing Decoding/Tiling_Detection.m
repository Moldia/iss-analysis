function Tiling_Detection(folder_image,image_prefix,image_suffix,...
    channel_max,x_size,y_size,show_tiled_position_YN,...
    low_resolution_full_size_image,create_CSV_file_YN,CSV_filename_prefix)

% re-tile images for detection (one hyb step)
%
% All images should be in the same folder.
% Image filename example: D:\colon\5plexdetection_c1_ORG.tif
%   folder_image = 'D:colon';
%   image_prefix = '5plexdetection';
%   image_suffix = '_ORG.tif';
%   channel_max = 6; (5 detection plus one for DAPI channel, if there is )
%
% Option1: create CSV file
%   Write csv file as CellProfiler input. 
%   The names in CSV file follow c1, c2,... as in the original image file.
%
% Option2: show tiled position
%   Draw yellow lines on a background image to show how the image is tiled.
%   Background image can be of any type, color/grayscale, and of low 
%   quality (preferred, save time for loading image). But the size has to 
%   be the same as original image.
%   Output image will be saved in the same folder as folder_image, with a
%   name of tiling_mxn.jpg (m: number of tiles in x dimension, n: number 
%   of tiles in y dimension)
%
% Tiling_Detection v2.1
% Apr 11, 2014, Xiaoyan

checkinput;
checkimages;

plot_x = zeros(ty+1,tx+1);
plot_y = zeros(ty+1,tx+1);
img_num = 0;
meta_pos = zeros(tx*ty,1);
tile_xpos = zeros(tx*ty,1);
tile_ypos = zeros(tx*ty,1);

% pad and tile images
%----------------------------
for c = 1:channel_max
    image_name = [folder_image image_prefix num2str(c) image_suffix];
    I = imread(image_name);
    img_num = img_num+1;
    if pad_x(c)
        zero_x = zeros(dim(c,2),pad_x(c));
        I = [I,zero_x];
    end
    if pad_y(c)
        zero_y = zeros(pad_y(c),dim(c,1)+pad_x(c));
        I = [I;zero_y];
    end
    % tile images
    %------------------------
    disp(['tiling image ' num2str(img_num) '..']);
    for i=1:ty
        for j=1:tx
            tile=I(y_size*(i-1)+1:y_size*i,x_size*(j-1)+1:x_size*j);
            out_dir = [folder_image 'Tiled_' image_prefix num2str(c) ...
                strtok(image_suffix,'.')];
            plot_x(i+1,j+1) = x_size*j;
            plot_y(i+1,j+1) = y_size*i;
            if exist(out_dir, 'dir')  
            else mkdir (out_dir);
            end
            tilenum = (i-1)*tx+j;
            imwrite(tile,[out_dir '\tile' num2str(tilenum) '.tif'],'tif');
            meta_pos(tilenum) = tilenum;
            tile_xpos(tilenum) = x_size*(j-1);
            tile_ypos(tilenum) = y_size*(i-1);
            path_name(tilenum,c) = {[out_dir '\']};
            file_name(tilenum,c) = {['tile' num2str(tilenum) '.tif']};
        end
    end
end

plot_x(1,2:end) = plot_x(2,2:end);
plot_y(2:end,1) = plot_y(2:end,2);


if show_tiled_position_YN
    showposition;
end
   
if create_CSV_file_YN
    writecsv;
end

% function to check input
    function checkinput
        if size(folder_image,2)
            folder_image =[folder_image '\'];
        end

        if show_tiled_position_YN 
            if exist(low_resolution_full_size_image,'file')==2
            else
                error('Could not find the image chosen to show tiled position.');
            end
        end
    end

% function to check input image size
    function checkimages
        dim = zeros(channel_max,2);
        disp('checking images..');
        for c = 1:channel_max
            image_name = [folder_image image_prefix num2str(c) image_suffix];
            if exist(image_name,'file')==2
                f = imfinfo(image_name);
            dim(c,:) =  [f.Width,f.Height];
            else
                error(['Could not find file ' image_name '.']);
            end 
        end

        max_x = max(dim(:,1));
        max_y = max(dim(:,2));
        tx = ceil(max_x/x_size);
        ty = ceil(max_y/y_size);

        new_x = x_size*tx;
        new_y = y_size*ty;
        pad_x = new_x - dim(:,1);
        pad_y = new_y - dim(:,2);
    end

% function write csv file
    function writecsv
        disp('writing CSV file..');
        fid = fopen([folder_image CSV_filename_prefix '.csv'], 'w');
        fprintf(fid,'Metadat_position,Tile_xPos,Tile_yPos,');
        for i = 1:channel_max
            fprintf(fid,'%s,%s,',['Image_PathName_c' num2str(i)],['Image_FileName_c' num2str(i)]);
        end
        fprintf(fid,'\n');
        for row = 1:tx*ty
            fprintf(fid,'%d,%d,%d,',meta_pos(row,:),tile_xpos(row),tile_ypos(row));
            for i = 1:channel_max
            fprintf(fid,'%s,%s,',path_name{row,i},file_name{row,i});
            end
            fprintf(fid,'\n');
        end 
        fclose(fid);
    end

% function to show tiled position
    function showposition
        disp('plotting to show tiled position..');
        Iref = imread(low_resolution_full_size_image);
        Iref = imresize(Iref,0.1);
        figure; imagesc(Iref); axis image; hold on;
        plot_x = floor(0.1*plot_x);
        plot_y = floor(0.1*plot_y);
        % vertical lines
        %--------------------
        for i = 1:tx+1
            x = plot_x(:,i);    
            y = plot_y(:,i);
            plot (x,y,'y');
        end
        % horizental lines
        %--------------------
        for i = 1:ty+1
            x = plot_x(i,:);    
            y = plot_y(i,:);
            plot (x,y,'y');
        end
        axis off;
        saveas(gcf, [folder_image 'tiling_' num2str(tx) 'x' num2str(ty)],'jpg');
    end
    
clear;
disp('Tiling finished.');
end

