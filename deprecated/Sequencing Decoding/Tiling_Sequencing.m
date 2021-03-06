function Tiling_Sequencing(folder_image,filename_base_prefix,in_subfolder_YN,...
    filename_channel_prefix,filename_suffix,base_max,channel_max,...
    x_size,y_size,standard_sequencing_protocol_YN,...
    show_tiled_position_YN,low_resolution_full_size_image,...
    separate_detection_step_YN,detection_image,...
    create_CSV_file_YN,CSV_filename_prefix,channel_order,varargin)

% Re-tile images for in situ sequencing image analysis.
% Run together with Sequencing_v2 script.
% 
% Image filename example1: D:\colon\base1_c1_ORG.tif
%   folder_image = 'D:colon';
%   filename_base_prefix = 'base';
%       in_subfolder_YN = 0;
%   filename_channel_prefix = '_c';
%   filename_suffix = '_ORG.tif';
%
% Image filename example2: D:\colon\base1\base1_c1_ORG.tif
%   folder_image = 'D:colon';
%   filename_base_prefix = 'base';
%       in_subfolder_YN = 1;
%   filename_channel_prefix = '_c';
%   filename_suffix = '_ORG.tif';
% 
% Choose a protocol:
%   Choice(1): Standard sequencing protocol 
%       The protocol used in most cases.
%   Choice(2): Separate detection step
%       The old protocol where detection of all blobs (general stain) 
%       is done in an independent step from sequencing.
%       The separate detection step image is given as additional input.
%       If there is also nuclei channel, reorder the files so that the base1
%       nuclei image is the one used in detection (so the blob-cell
%       relation will be always based on the detection step.
%
% Option1: Create CSV file
%   Write csv file as CellProfiler input. 
%   When specifying the channel order, use standarized names, 
%       eg. channel_order = {'Nuclei' 'General_stain' 'T' 'G' 'C' 'A'}; 
%       especially for Nuclei and General_stain.
%   The same channel order should be strictly followed in every imaging
%   cycle. If not, manually change the original image filenames.
%   There is no need to have every one of them (even applicable for Nuclei
%   and General_stain).
%   If using the alternative protocol (separate detection), do not include
%   'General_stain' in channel_order. It will be added automatically.
%
% Option2: Show tiled position
%   Draw yellow lines on a background image to show how the image is tiled.
%   Background image can be of any type, color/grayscale, and of low 
%   quality (preferred, save time for loading image). But the size has to 
%   be the same as original image.
%   Output image will be saved in the same folder as folder_image, with a
%   name of tiling_mxn.jpg (m: number of tiles in x dimension, n: number 
%   of tiles in y dimension)
%
%
% Tiling_Sequencing v2.1.2
% Aug 6, 2014, Xiaoyan

checkinput;   
checkimages;

% predefine the size of arrays
%----------------------------
nonuclei = 0;
nogeneral = 0;
plot_x = zeros(ty+1,tx+1);
plot_y = zeros(ty+1,tx+1);
img_num = 0;
meta_pos = zeros(tx*ty*base_max,1);
tile_xpos = zeros(tx*ty*base_max,1);
tile_ypos = zeros(tx*ty*base_max,1);
hyb_step = cell(tx*ty*base_max,1);
path_name = cell(tx*ty*base_max,c);
file_name = cell(tx*ty*base_max,c);
general_blob = cell(tx*ty*base_max,2);
nuclei_file = cell(tx*ty*base_max,2);

% pad and tile images
%----------------------------
for b = 1:base_max
    for c = 1:channel_max
        if in_subfolder_YN
            image_name = [folder_image filename_base_prefix ...
                num2str(b) '\' filename_base_prefix num2str(b) ...
                filename_channel_prefix num2str(c) filename_suffix];
        else
            image_name = [folder_image filename_base_prefix  num2str(b) ...
                filename_channel_prefix num2str(c) filename_suffix];
        end
        I = imread(image_name);
        img_num = img_num+1;
        padx = pad_x((b-1)*channel_max+c+separate_detection_step_YN);
        pady = pad_y((b-1)*channel_max+c+separate_detection_step_YN);
        if padx
            dimy = dim((b-1)*channel_max+c+separate_detection_step_YN,2);
            zero_x = zeros(dimy,padx);
            I = [I,zero_x];
        end
        if pady
            dimx = dim((b-1)*channel_max+c+separate_detection_step_YN,1);
            zero_y = zeros(pady,dimx+padx);
            I = [I;zero_y];
        end
        % tile images
        %------------------------
        disp(['tiling image ' num2str(img_num) '..']);
        for i=1:ty
            for j=1:tx
                tilenum = (i-1)*tx+j;
                if skip_image_tiling_YN
                    if in_subfolder_YN
                        out_dir = [folder_image filename_base_prefix ...
                        num2str(b)];
                        path_name((tilenum-1)*base_max+b,c) = {[out_dir '\']};
                      else
                        path_name((tilenum-1)*base_max+b,c) = {folder_image};
                    end
                    file_name((tilenum-1)*base_max+b,c) = {[filename_base_prefix num2str(b) ...
                        filename_channel_prefix num2str(c) filename_suffix]};
                else
                    tile=I(y_size*(i-1)+1:y_size*i,x_size*(j-1)+1:x_size*j);
                    out_dir = [folder_image 'Tiled_' filename_base_prefix ...
                        num2str(b) filename_channel_prefix num2str(c) ...
                        strtok(filename_suffix,'.')];

                    if exist(out_dir, 'dir')  
                        else mkdir (out_dir);
                    end
                    plot_x(i+1,j+1) = x_size*j;
                    plot_y(i+1,j+1) = y_size*i;

                    imwrite(tile,[out_dir '\tile' num2str(tilenum) '.tif'],'tif');
                    path_name((tilenum-1)*base_max+b,c) = {[out_dir '\']};
                    file_name((tilenum-1)*base_max+b,c) = {['tile' num2str(tilenum) '.tif']};
                end
                if c ==1
                    meta_pos((tilenum-1)*base_max+b) = tilenum;
                    tile_xpos((tilenum-1)*base_max+b) = x_size*(j-1);
                    tile_ypos((tilenum-1)*base_max+b) = y_size*(i-1);
                    hyb_step((tilenum-1)*base_max+b) = {['hyb' num2str(b)]};
                end
            end
        end
    end
end
plot_x(1,2:end) = plot_x(2,2:end);
plot_y(2:end,1) = plot_y(2:end,2);


if standard_sequencing_protocol_YN
    generalblob;
elseif separate_detection_step_YN
    detection;
end

nuclei;

if show_tiled_position_YN
    showposition
end

if create_CSV_file_YN
    writecsv;
end


% function to check input
    function checkinput
        if size(folder_image,2)
            if folder_image(end)=='\'
            else
                folder_image =[folder_image '\'];
            end
        end

        if show_tiled_position_YN 
            if exist(low_resolution_full_size_image,'file')==2
            else
                error('Could not find the image chosen to show tiled position.');
            end
        end

        if separate_detection_step_YN == 0 || separate_detection_step_YN ==1
        else
            error('separate_detection_step_YN should be set either to 0 or to 1.');
        end

        if standard_sequencing_protocol_YN
        elseif separate_detection_step_YN
        else
            error(['Currently only supports standard protocol and '...
                'separate detection step two types of sequencing protocols. '...
                'Choose at least one of them.']);
        end
        
        if create_CSV_file_YN
            if length(channel_order) == channel_max
            else
                error(['Maximum channel number and the number of elements '...
                    'in channel_order does not match.']);
            end
        end
        
        switch length(varargin)
            case 0
                skip_image_tiling_YN = 0;
            case 1
                skip_image_tiling_YN = varargin{1};
        end
    end

% function to check input image size
    function checkimages
        dim = zeros(base_max*channel_max+separate_detection_step_YN,2);
        disp('checking images..');
        if separate_detection_step_YN
            if exist(detection_image,'file')==2
                f = imfinfo(detection_image);
                dim(1,:) = [f.Width,f.Height];
            else
                error('Could not find the blob detection image.');
            end
        end
        for b = 1:base_max
            for c = 1:channel_max
                if in_subfolder_YN
                    image_name = [folder_image filename_base_prefix ...
                        num2str(b) '\' filename_base_prefix num2str(b) ...
                        filename_channel_prefix num2str(c) filename_suffix];
                else
                    image_name = [folder_image filename_base_prefix  num2str(b) ...
                        filename_channel_prefix num2str(c) filename_suffix];
                end
                if exist(image_name,'file')==2
                    f = imfinfo(image_name);
                    dim((b-1)*channel_max+c+separate_detection_step_YN,:) = [f.Width,f.Height];
                else
                    error(['Could not find file ' image_name '.']);
                end 
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

% function to fix general blob
    function generalblob
        bind = find(ismember(channel_order,...
            {'General_stain','general_stain','General_blob','general_blob',...
            'General stain','general stain','General blob','general blob',...
            'AF750','Cy7'}));
        if length(bind) == 1
            channel_order(bind) = {'Spec_blob'};
            for i = 1:tx*ty*base_max
                switch mod(i,base_max)
                    case 1
                        general_blob(i,1) = path_name(i,bind);
                        general_blob(i,2) = file_name(i,bind);
                    case 0
                        general_blob(i,1) = general_blob(i-base_max+1,1);
                        general_blob(i,2) = general_blob(i-base_max+1,2);
                    otherwise
                        general_blob(i,1) = general_blob(i-mod(i,base_max)+1,1);
                        general_blob(i,2) = general_blob(i-mod(i,base_max)+1,2);
                end
            end
        elseif isempty(bind)
            nogeneral = 1;
        else
            error('Specify only one General_stain in channel_order.');
        end
    end

% function to tile and fix general blob in separate detection protocol
    function detection
        I = imread(detection_image);
            img_num = img_num+1;
            padx = pad_x(1);
            pady = pad_y(1);
            if padx
                dimy = dim(1,2);
                zero_x = zeros(dimy,padx);
                I = [I,zero_x];
            end
            if pady
                dimx = dim(1,1);
                zero_y = zeros(pady,dimx+padx);
                I = [I;zero_y];
            end

        disp(['tiling image ' num2str(img_num) '..']);
        for i=1:ty
            for j=1:tx
                tilenum = (i-1)*tx+j;
                if skip_image_tiling_YN
                    spli = strsplit(detection_image,'\');
                    out_dir = strjoin(spli(1:end-1),'\');
                    general_blob((tilenum-1)*base_max+1:tilenum*base_max,1) = {[out_dir '\']};
                    general_blob((tilenum-1)*base_max+1:tilenum*base_max,2)...
                        = spli(end);
                else
                    tile=I(y_size*(i-1)+1:y_size*i,x_size*(j-1)+1:x_size*j);
                    out_dir = [folder_image 'Tiled_detection'];
                    if exist(out_dir, 'dir')  
                    else mkdir (out_dir);
                    end
                    imwrite(tile,[out_dir '\tile' num2str(tilenum) '.tif'],'tif');
                    general_blob((tilenum-1)*base_max+1:tilenum*base_max,1) = {[out_dir '\']};
                    general_blob((tilenum-1)*base_max+1:tilenum*base_max,2)...
                        = {['tile' num2str(tilenum) '.tif']};
                end
            end
        end
    end
    
% function to fix nuclei
    function nuclei
        nind = find(ismember(channel_order,{'Nuclei','nuclei','DAPI','Hoechst'}));
        if length(nind) == 1
            channel_order(nind) = [];
            for i = 1:tx*ty*base_max
                switch mod(i,base_max)
                    case 1
                        nuclei_file(i,1) = path_name(i,nind);
                        nuclei_file(i,2) = file_name(i,nind);
                    case 0
                        nuclei_file(i,1) = nuclei_file(i-base_max+1,1);
                        nuclei_file(i,2) = nuclei_file(i-base_max+1,2);
                    otherwise
                        nuclei_file(i,1) = nuclei_file(i-mod(i,base_max)+1,1);
                        nuclei_file(i,2) = nuclei_file(i-mod(i,base_max)+1,2);
                end
            end
            path_name(:,nind) = [];
            file_name(:,nind) = [];
        elseif isempty(nind)
            nonuclei = 1;
        else
            error('Specify only one Nuclei in channel_order.');
        end
    end
 
% function write csv file
    function writecsv
        disp('writing CSV file..');
        % replace all backslash
%         path_name = strrep(path_name,'\','/');

        fid = fopen([folder_image CSV_filename_prefix '.csv'], 'w');
       
        % uniform header
        fprintf(fid,'%s,%s,%s,%s','Metadata_position',...
                'Tile_xPos','Tile_yPos','Hyb_step');
        if nogeneral
            % nuclei header
            if nonuclei
            else
                fprintf(fid,',Image_PathName_Nuclei,Image_FileName_Nuclei');
            end
            % channel header
            for i = 1:length(channel_order)
                    fprintf(fid,',%s,%s',['Image_PathName_' channel_order{i}],...
                        ['Image_FileName_' channel_order{i}]);
            end
            fprintf(fid,'\n');
            
            % metadata, tile position, hyb step and nuclei(if any) columns
            for row = 1:tx*ty*base_max
                if nonuclei
                    % uniform info
                    fprintf(fid,'%d,%d,%d,%s',meta_pos(row),...
                        tile_xpos(row),tile_ypos(row),hyb_step{row});
                else
                    % replace backslash
%                     nuclei_file(:,1) = strrep(nuclei_file(:,1),'\','/');
                    % uniform info + nuclei
                    fprintf(fid,'%d,%d,%d,%s,%s,%s',meta_pos(row),...
                        tile_xpos(row),tile_ypos(row),hyb_step{row},...
                        nuclei_file{row,1},nuclei_file{row,2});
                end
                % channel info
                for i = 1:length(channel_order)
                    fprintf(fid,',%s,%s',path_name{row,i},file_name{row,i});
                end
                fprintf(fid,'\n');
            end
        else
            % general blob header
            fprintf(fid,',%s,%s',...
                'Image_PathName_General_blob','Image_FileName_General_blob');
            if nonuclei
            else
                 fprintf(fid,',Image_PathName_Nuclei,Image_FileName_Nuclei');
            end
            % channel header
            for i = 1:length(channel_order)
                fprintf(fid,',%s,%s',['Image_PathName_' channel_order{i}],['Image_FileName_' channel_order{i}]);
            end
            fprintf(fid,'\n');
            for row = 1:tx*ty*base_max
                if nonuclei
%                     general_blob(:,1) = strrep(general_blob(:,1),'\','/');
                    fprintf(fid,'%d,%d,%d,%s,%s,%s,',meta_pos(row),...
                        tile_xpos(row),tile_ypos(row),hyb_step{row},...
                        general_blob{row,1},general_blob{row,2});
                else
%                     nuclei_file(:,1) = strrep(nuclei_file(:,1),'\','/');
%                     general_blob(:,1) = strrep(general_blob(:,1),'\','/');
                    fprintf(fid,'%d,%d,%d,%s,%s,%s,%s,%s',meta_pos(row),...
                        tile_xpos(row),tile_ypos(row),hyb_step{row},...
                        general_blob{row,1},general_blob{row,2},...
                        nuclei_file{row,1},nuclei_file{row,2});
                end
                for i = 1:length(channel_order)
                    fprintf(fid,',%s,%s',path_name{row,i},file_name{row,i});
                end
                fprintf(fid,'\n');
            end
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

