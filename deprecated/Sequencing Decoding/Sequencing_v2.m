% Integrated analysis for sequencing.
% For a specific module, if you don't know how to set the parameteres,
% eg. Tiling_Sequencing, type "help Tiling_Sequencing" in the command
% window.
% Ideally keep one file for one experiment
% Five major functions: 1)Tiling, 2)Decoding, 3)Prelimiary analysis before QT,
% 4)Thresholding, 5)Plotting (global or tile)
% Carefully set the paramets.
% all variables ending with _YN are yes or no (1 or 0) questions
%
% Sequencing v2.4.1
% Xiaoyan, 2015-5-3


clear, clc; close all; drawnow;

% Choose functions to run
run_Tiling_YN = 0;
run_Decode_YN = 0;
run_Analysis_YN = 0;
run_Threshold_YN = 0;
run_Plotting_Global_YN = 0;
run_Plotting_Tile_YN = 0;
%================================================
% set parameters
%----------------------------
% Tiling_Sequencing
    folder_image = 'D:\experiment\colon'; % preferably full path name
    filename_base_prefix = 'base';  % keep single quote marks
        in_subfolder_YN = 0;
    filename_channel_prefix = '_c';
    filename_suffix = '_ORG.tif';
    base_max = 4;       channel_max = 6;
    x_size = 3000;      y_size = 3000;
    % choose a protocol
    standard_sequencing_protocol_YN = 1;
    separate_detection_step_YN = 0;
        detection_image = ''; % preferably full path
    % options
    create_CSV_file_YN = 1; 
        CSV_filename_prefix = 'Tiled';
        channel_order = {'Nuclei' 'T' 'G' 'C' 'A' 'General_stain'};
    show_tiled_position_YN = 0; 
        low_resolution_full_size_image = 'C:\experiment\colon\background.tif';
    skip_image_tiling_YN = 0; % if you only want the csv file or the image showing tiled position (of course the tiled images should already be there)
%----------------------------
% Decoding_Sequencing
% if you did Tiling using this script, first move the script file to the folder that contains CP results
    input_file = 'blobs.csv';
    % don't change this unless your input file has some weird form
    General_Alignment_ACGT_column_number = [4,5,6,7,8,9];    % use 0 if any of them is MISSING in the file
    XYPosition_ParentCell_column_number = [10,11,12];
    
    num_hybs = 4;  % if 5th empty cycle exists, num_hybs=5, no upper limit
        cycle5_empty_YN = 0;
    taglist = ID_empty;
        use_old_style_taglist_YN = 0;  % the digit column will be simply ignored
        grouped_YN = 0;
    calculate_global_position_YN = 1;
        retiled_by_Python_or_MATLAB_YN = 1;
            csv_file_contain_tile_position = 'CSVmodified.csv'; %full path if not in the current directory
        retiled_by_ZEN_YN = 0;  % may have a few pixels shift
            tile_x_dimension_max = 2000;    x_tile_number = 10; % tiles along x-axis, in column
            tile_y_dimension_max = 2000;    y_tile_number = 10; % tiles along y-axis, in row
    output_directory_decode = 'Decoding';   
    output_filename_decode_prefix = 'beforeQT';
    % options
    check_parent_cell_YN = 0;       
    check_alignment_YN = 0;
        alignment_min_threshold = 1.8;
    remove_T_in_analysis = 0;
    abnormal_sequencing_YN = 0;
        sequencing_order = '1234';  % keep the quote marks, same length as (num_hybs - cycle5_empty_YN)
%----------------------------
% Analysis_Sequencing_beforeQT
    narrow_down_quality_range_YN = 0;
        lower_limit = 0.3;        upper_limit = 0.6;
    table_gene_counts_at_different_thresholds_YN = 0;
        threshold = [0.1, 0.3, 0.4, 0.8];
    plot_general_quality_items_YN = 1;
    plot_seq_spec_channel_histogram_YN = 0;
    guess_closest_expected_reads_YN = 0;
    plot_quality_vs_general_stain_YN = 0;
    plot_spectrum_YN = 0;
%----------------------------
% Threshold_Sequencing
    quality_threshold = 0.4;        general_strain_threshold = 0.01;
    % saved in the same folder as beforeQT, default name: QT_quality_generalstain
    output_filename_afterQT_prefix = ['QT_' num2str(quality_threshold) '_' num2str(general_strain_threshold)];
    % options
    remove_homomer_YN = 0;
        based_on_empty_cycle_5_YN = 0; % of course you should have settings correct in beforeQT decoding
            remove_only_above_quality_threshold_YN = 0; 
%----------------------------
% Plotting_global_Sequencing
    background_image = 'aligned_HEstain_10%.jpg'; 
        reduced_size = 0.1;
        I_want_to_plot_on_white_backgound = 0;
    taglist_plot = taglist; % default: use the same one as in Decoding
        use_old_style_YN = 0;  % digit column will be ignored
        use_default_symbol_list_YN = 0;
    symbol_size = 8; % default: 6
    output_directory_plot = 'Plotting';     
    output_filename_plot_prefix = 'AllPlottedOnHE';
    % options
    exclude_NNNN_YN = 0;
    plot_reads_beforeQT_YN = 0;
    plot_based_on_group_YN = 0;
    plot_base1_general_stain = 0;
%----------------------------
% Plotting_tile_sequencing
    single_tile_YN = 0;
        tile_background_image = 'D:\experiment\colon\singleTile_HEstain.jpg'; % full path
    multiple_tiles_YN = 1;
        background_image_tile_prefix = 'D:\experiment\colon\HEstain\HE_Tile'; %include the full path
        background_image_tile_suffix = '.jpg';
            continuous_tile_id_YN = 1; % from tile id 1 to max, all will be plotted
                max_tile_number = 20;
            subset_of_tiles_YN = 0; % only the defined tiles will be plotted
                tile_id_to_plot = [1 3]; 
    taglist_plot_tile = taglist; % default: use the same one as in Decoding
        tile_use_old_style_taglist_YN = 0;
        tile_use_default_symbol_list = 0;
    symbol_size_tile = 8; % default: 6
    output_directory_plot_tile = 'Plotting';     
    output_filename_plot_tile_prefix = 'AllPlottedOnTile';
    % options
    plot_tile_exclude_NNNN_YN = 0;
    plot_tile_reads_beforeQT_YN = 0;
    plot_tile_based_on_group_YN = 0;
    plot_tile_base1_general_stain = 0;

 
%================================================
if run_Tiling_YN || run_Decode_YN || run_Analysis_YN || run_Threshold_YN || run_Plotting_Global_YN || run_Plotting_Tile_YN
else
    error('Choose at least one function.');
end

if run_Tiling_YN
    Tiling_Sequencing(folder_image,filename_base_prefix,in_subfolder_YN,...
        filename_channel_prefix,filename_suffix,base_max,channel_max,...
        x_size,y_size,standard_sequencing_protocol_YN,...
        show_tiled_position_YN,low_resolution_full_size_image,...
        separate_detection_step_YN,detection_image,...
        create_CSV_file_YN,CSV_filename_prefix,channel_order,...
        skip_image_tiling_YN);
end
if run_Decode_YN
    Decoding_Sequencing(input_file,num_hybs,cycle5_empty_YN,...
        taglist,use_old_style_taglist_YN,grouped_YN,...
        calculate_global_position_YN,retiled_by_Python_or_MATLAB_YN,...
        csv_file_contain_tile_position,retiled_by_ZEN_YN,...
        tile_x_dimension_max,tile_y_dimension_max,x_tile_number,y_tile_number,...
        output_directory_decode,output_filename_decode_prefix,...
        check_parent_cell_YN,check_alignment_YN,alignment_min_threshold,...
        remove_T_in_analysis,abnormal_sequencing_YN,sequencing_order,...
        General_Alignment_ACGT_column_number,XYPosition_ParentCell_column_number);
end
if run_Analysis_YN
    Analysis_Sequencing_beforeQT(output_directory_decode,...
        output_filename_decode_prefix,...
        narrow_down_quality_range_YN,lower_limit,upper_limit,...
        table_gene_counts_at_different_thresholds_YN,threshold,...
        plot_general_quality_items_YN,plot_seq_spec_channel_histogram_YN,...
        guess_closest_expected_reads_YN,plot_quality_vs_general_stain_YN,...
        plot_spectrum_YN);
end
if run_Threshold_YN
    Threshold_Sequencing(output_directory_decode,...
        output_filename_decode_prefix,output_filename_afterQT_prefix,...
        quality_threshold,general_strain_threshold,remove_homomer_YN,...
        based_on_empty_cycle_5_YN,remove_only_above_quality_threshold_YN);
end 
if run_Plotting_Global_YN
    Plotting_global_Sequencing(output_directory_decode,...
        output_filename_afterQT_prefix,background_image,reduced_size,...
        num_hybs,cycle5_empty_YN,abnormal_sequencing_YN,...
        sequencing_order,taglist_plot,use_old_style_YN,...
        use_default_symbol_list_YN,symbol_size,...
        plot_reads_beforeQT_YN,output_directory_plot,...
        output_filename_plot_prefix,exclude_NNNN_YN,...
        plot_based_on_group_YN,...
        plot_base1_general_stain,I_want_to_plot_on_white_backgound)
end
if run_Plotting_Tile_YN
    Plotting_tile_Sequencing(output_directory_decode,...
        output_filename_afterQT_prefix,...
        single_tile_YN,tile_background_image,...
        multiple_tiles_YN,background_image_tile_prefix,background_image_tile_suffix,...
        continuous_tile_id_YN,max_tile_number,subset_of_tiles_YN,tile_id_to_plot,...
        taglist_plot_tile,tile_use_old_style_taglist_YN,symbol_size_tile,...
        output_directory_plot_tile,output_filename_plot_tile_prefix,...
        plot_tile_exclude_NNNN_YN,plot_tile_reads_beforeQT_YN,...
        plot_tile_based_on_group_YN,tile_use_default_symbol_list,...
        plot_tile_base1_general_stain)
end
%----------------------------
% fig = findobj;
% nf = max(fig(fig==fix(fig)));
% while nf >= 1
%     figure(nf);
%     nf = nf-1;
% end

clear;   
    