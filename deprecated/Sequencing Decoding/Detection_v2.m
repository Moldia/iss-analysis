%********************************************************************
% Tools for detection
% Sep 21, 2014, Xiaoyan
%********************************************************************
clear, clc;

% Choose functions to Derun
run_Tiling_YN = 1;

%================================================
% set parameters
%----------------------------
% Tiling_detection
if run_Tiling_YN
    folder_image= 'D:\Degree_project\Imaging\140213_canceroid\reverse_GCAT_b1-5\ZENexport'; % preferably full path name
    image_prefix = 'base2_c';  % keep single quote marks
    image_suffix = '_ORG.tif';
    channel_max = 6;
    x_size = 3000;      y_size = 3000;
    % options
    create_CSV_file_YN = 1;
        CSV_filename_prefix = 'CPinput_3';
     show_tiled_position_YN = 1; 
        low_resolution_full_size_image = 'D:\Degree_project\Imaging\140213_canceroid\reverse_GCAT_b1-5\base1_display_c2.jpg';
end
%----------------------------

%================================================
% if run_Tiling_YN || run_Decode_YN || run_Analysis_YN || run_Threshold_YN || run_Plotting_Global_YN || run_Plotting_Tile_YN
% else
%     error('Choose at least one function.');
% end

if run_Tiling_YN
    Tiling_Detection(folder_image,image_prefix,image_suffix,...
        channel_max,x_size,y_size,show_tiled_position_YN,...
        low_resolution_full_size_image,create_CSV_file_YN,CSV_filename_prefix);
end


clear;   
    