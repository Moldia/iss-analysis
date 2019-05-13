function  Decoding_Sequencing(input_file,num_hybs,cycle5_empty_YN,...
    taglist,use_old_style_taglist_YN,grouped_YN,...
    calculate_global_position_YN,retiled_by_Python_or_MATLAB_YN,...
    csv_file_contain_tile_position,retiled_by_ZEN_YN,...
    tile_x_dimension_max,tile_y_dimension_max,x_tile_number,y_tile_number,...
    output_directory_decode,output_filename_decode_prefix,...
    check_parent_cell_YN,check_alignment_YN,alignment_min_threshold,...
    remove_T_in_analysis,abnormal_sequencing_YN,sequencing_order,varargin)

% Decode the CellProfile (CP) result to reads
% The columns in the CP output must be 1.ImageNumber 2.ObejctNumber
% 3. MetadataPosition 4.Intensity_general_stain 5.Intensity_ImageMath
% 6-9.Intensity_A,C,G,T 10-11.x,y_postion 12.Parent_Cell (12 can be
% missing)
%
% Taglist settings:
%   Automatically converts barcode letters to digits.
%   Choice(1): Use old style taglist. 
%       Eg. 'ATTA HES1', 4114,'r.';
%       The taglist used before. Make sure the barcodes are correct. 
%       The digit column will be sipmly ignored.
%   Choice(2): Group the tags and give group tag.
%       Eg. 'ATCG STK15 mutation','go';
%       The same group tag can be shared by several barcodes.
%
% Calculate global position
%   Choice(1): Re-tiled by Python or MATLAB.
%       Specify the CSV file used as CellProfiler input. It contains the
%       x- and y- starting position of tiles.
%   Choice(2): Re-tiled by ZEN.
%       Manually find the biggest x- and y- size of tiled images. Also set
%       how many tiles are in x- and y- dimension.
%       The global position will be calculated by adding up the tile sizes.
%
% Ouput directoy: a sub-folder will be created within the current working
% directory.
%
% Option1: Check parent cell.
%   All blobs without a parent cell will be discarded.
%   Set to 1 if the blob-cell relation is fairly good in CellProfiler.
%
% Option2: Check alignment.
%   The recommended alignemnt score threshold is 1.8-2
%   Only blobs with good alignemnt in all hybridization steps will be
%   saved. Others will be discarded.
%   Not always necessary, since the general stain for each hybridization
%   step will be checked later in thresholding.
%
% Option3: Remove T in analysis.
%   The T channel will be directly excluded in decoding. 
%   But an additional read guess regarding T will be done and give a ui
%   table.
%   For now, these reads are still considered as unexpected.
%
% Option4: Abnormal sequencing.
%   If you use some special sequencing order, 
%   or if something has gone wrong in CellProfiler analysis...
%   The length of sequencing order has to match with number of hybs (except
%   when there is 5th empty cycle). The number at each position should be
%   bettween 0 to 4 (base position in barcode).
%   Example1: reverse sequencing.
%       sequencing_order = '4321';
%   Example2: sequenced base3 twice and you want reads from both
%       num_hybs = 5; sequencing_order = '12334';
%   Example3: sequenced base3 twice and you want read only from the second 
%       num_hybs = 5; sequencing_order = '12034';
%   Example4: base2 has very poor alignment and you want exclude it
%       num_hybs = 4; sequencing_order = '1034';
%   
%
% Decoding_Sequencing v2.4.2
% Xiaoyan, 2015-5-3

%% initiate
drawnow;
disp('Initiating Decoding_Sequencing.');
total_t1 = clock;

%% import and structure data, reform taglist
checkinputandimport;
if remove_T_in_analysis
    seqdata(:,9) = 1E-6;
end
[taglist,expected_list,expected_digits,~,exp_tags,exp_groups] =  ...
    maketag(taglist,num_hybs,cycle5_empty_YN,use_old_style_taglist_YN,...
    grouped_YN,abnormal_sequencing_YN,sequencing_order);

%% preallocate arrays
seq_res = zeros(num_blobs,num_hybs); 
channel_strength_max = zeros(num_blobs,num_hybs);
general_strength_sum = zeros(num_blobs,num_hybs);
num_channels = 4;
channel_strength_original = zeros(num_hybs,num_channels,num_blobs);
seq_quality = zeros(num_blobs,num_hybs);
general_strength = zeros(num_blobs,num_hybs);
alignment_score = zeros(num_blobs,num_hybs);
tile_ID = zeros(num_blobs,1);
tile_cell_ID = zeros(num_blobs,1);
cell_ID = zeros(num_blobs,1);
tile_blob_ID = zeros(num_blobs,1);
blob_ID = zeros(num_blobs,1);
tile_x_pos = zeros(num_blobs,1);
tile_y_pos = zeros(num_blobs,1);

if calculate_global_position_YN
    global_x_pos = zeros(num_blobs,1);
    global_y_pos = zeros(num_blobs,1);
end


%% main function: extract data
start_blob_ID=0;
start_cell_ID=0;
fprintf('number of tiles\t number of blobs\n'); 
fprintf('\t%6d\t\t%6d\n',num_tiles,num_blobs);
fprintf('# tile\tnumber of blobs\n');

for t=1:num_tiles
    temp_tile_seq_data = seqdata(seqdata(:,3) == t,:);
    temp_tile_num_blobs = max(temp_tile_seq_data(:,2));
    temp_tile_num_cells = max(temp_tile_seq_data(:,12));
    
    %------------------------
    if calculate_global_position_YN && retiled_by_Python_or_MATLAB_YN
        temp_tile_x_start = mean(tile_start_pos(tile_start_pos(:,1) == t,2));
        temp_tile_y_start = mean(tile_start_pos(tile_start_pos(:,1) == t,3));
    elseif calculate_global_position_YN && retiled_by_ZEN_YN
        if mod(t,x_tile_number) == 0
            temp_tile_y_start = (floor(t/x_tile_number)-1)*tile_y_dimension_max;
            temp_tile_x_start = (x_tile_number-1)*tile_x_dimension_max;
        else 
            temp_tile_y_start = floor(t/x_tile_number)*tile_y_dimension_max;
            temp_tile_x_start = (mod(t,x_tile_number)-1)*tile_x_dimension_max;
        end
    end
    fprintf('%6d\t\t%6d\n',t,temp_tile_num_blobs);

    %------------------------
    if temp_tile_num_blobs
        [~,order] = sort(temp_tile_seq_data(:,2));
                
        % extract info about blobs within the tile
        temp_tile_intensity = temp_tile_seq_data(order,6:6+num_channels-1);
        temp_tile_general = temp_tile_seq_data(order,4);
        temp_tile_align = temp_tile_seq_data(order,5);
        temp_tile_x_pos = temp_tile_seq_data(order,10);
        temp_tile_y_pos = temp_tile_seq_data(order,11);
        temp_tile_cell_id = temp_tile_seq_data(order,12);
        
        % max intensities
        [maxint,maxidx] = max(temp_tile_intensity,[],2);
        channel_strength_max(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(maxint,num_hybs,temp_tile_num_blobs))';
        
        % seq result matrix
        seq_res(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(maxidx,num_hybs,temp_tile_num_blobs))';
        
        % blob info
        tile_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = t;
        tile_blob_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            (1:temp_tile_num_blobs)';
        blob_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            (1:temp_tile_num_blobs)'+repmat(start_blob_ID,temp_tile_num_blobs,1);
        alignment_score(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(temp_tile_align,num_hybs,temp_tile_num_blobs))';

        
        temp_list = 1:temp_tile_num_blobs*num_hybs;
        if num_hybs>1
            temp_list = mod(temp_list,num_hybs)==1;  
        end
        
        tile_cell_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_cell_id(temp_list);
        cell_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_cell_id(temp_list)+(temp_tile_cell_id(temp_list)~=0)*start_cell_ID;
        tile_x_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_x_pos(temp_list);
        tile_y_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_y_pos(temp_list);
        % global coordinates
        if calculate_global_position_YN
            global_x_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
                repmat(temp_tile_x_start,temp_tile_num_blobs,1) + temp_tile_x_pos(temp_list);
            global_y_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
                repmat(temp_tile_y_start,temp_tile_num_blobs,1) + temp_tile_y_pos(temp_list);
        end
       
        % original intensities
        for c = 1:num_channels
            channel_strength_original(:,c,start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
                reshape(temp_tile_intensity(:,c),num_hybs,1,temp_tile_num_blobs);
        end
        general_strength(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
           (reshape(temp_tile_general,num_hybs,temp_tile_num_blobs))';
        general_strength_sum(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            reshape(sum(temp_tile_intensity,2),num_hybs,temp_tile_num_blobs)';
        
        start_blob_ID = start_blob_ID+temp_tile_num_blobs;
        start_cell_ID = start_cell_ID+temp_tile_num_cells;
    end

end
clear -regexp ^temp
clear t c

%% bases into reads
[seq_num,inuni] = combinebase(abnormal_sequencing_YN,num_hybs,cycle5_empty_YN,...
    num_blobs,seq_res,sequencing_order);

%% sequencing quality
seq_quality = channel_strength_max./general_strength_sum;
qualityitems;
qualitycheck;

%% all reads before QT
allbt = seq_num(ffbt); 

if abnormal_sequencing_YN
    seq_res_allbt = seq_res(ffbt,inuni);
    seq_quality_allbt = seq_quality(ffbt,inuni);
    general_strength_allbt = general_strength(ffbt,inuni);
    alignment_score_allbt = alignment_score(ffbt,inuni);
    channel_strength_max_allbt = channel_strength_max(ffbt,inuni);
    channel_strength_original_allbt = channel_strength_original(inuni,:,ffbt);   
    if cycle5_empty_YN
        seq_res_b5_bt = seq_res_b5(ffbt);
        seq_quality_b5_bt = seq_quality_b5(ffbt);
        general_strength_b5_bt = general_strength_b5(ffbt);
        channel_max_b5_bt = channel_max_b5(ffbt);
    end
else
    if cycle5_empty_YN
        seq_res_allbt = seq_res(ffbt,1:4);
        seq_quality_allbt = seq_quality(ffbt,1:4);
        general_strength_allbt = general_strength(ffbt,1:4);
        alignment_score_allbt = alignment_score(ffbt,1:4);
        channel_strength_max_allbt = channel_strength_max(ffbt,1:4);
        channel_strength_original_allbt = channel_strength_original(1:4,:,ffbt);
        
        seq_res_b5_bt = seq_res_b5(ffbt);
        seq_quality_b5_bt = seq_quality_b5(ffbt);
        general_strength_b5_bt = general_strength_b5(ffbt);
        channel_max_b5_bt = channel_max_b5(ffbt);
    else
        seq_res_allbt = seq_res(ffbt,:);
        seq_quality_allbt = seq_quality(ffbt,:);
        general_strength_allbt = general_strength(ffbt,:);
        alignment_score_allbt = alignment_score(ffbt,:);
        channel_strength_max_allbt = channel_strength_max(ffbt,:);
        channel_strength_original_allbt = channel_strength_original(:,:,ffbt);
    end
end
letters_allbt = num2letters_2(seq_res_allbt);

tile_allbt = tile_ID(ffbt);
tile_cell_allbt = tile_cell_ID(ffbt);
cell_allbt = cell_ID(ffbt);
tile_blob_allbt = tile_blob_ID(ffbt);
blob_allbt = blob_ID(ffbt);
tile_x_allbt = tile_x_pos(ffbt);
tile_y_allbt = tile_y_pos(ffbt);
if calculate_global_position_YN
    global_x_allbt = global_x_pos(ffbt);
    global_y_allbt = global_y_pos(ffbt);
end

seq_quality_min_allbt = seq_quality_min(ffbt,:);
general_strength_min_allbt = general_strength_min(ffbt,:);
alignment_score_min_allbt = alignment_score_min(ffbt,:);
seq_quality_min_ind_allbt = seq_quality_min_ind(ffbt,:);
general_strength_min_ind_allbt = general_strength_min_ind(ffbt,:);
alignment_score_min_ind_allbt = alignment_score_min_ind(ffbt,:);

%% assign gene name and group tag to each read
[name_tag_allbt,group_tag_allbt] = ...
    genename(grouped_YN,allbt,expected_list,exp_tags,exp_groups);

if remove_T_in_analysis
    [Tdata,Tcolumnn,Trown] = ...
        GuessT(expected_digits(:,2:end),allbt,name_tag_allbt,taglist(:,1));
end

%% group reads
[uni_seq_res_allbt,count_code_allbt,uni_letters_allbt,uni_reads_name_tag_allbt,uni_reads_group_tag_allbt,catagory_read_allbt] = ...
    groupreads(allbt,seq_res_allbt,letters_allbt,expected_list,exp_tags,exp_groups);

% find unexpected homomer reads
homomer_index_allbt = findHomomer(expected_list,allbt);
catagory_read_allbt(homomer_index_allbt) = -1;

%% output figures and files
disp('busy drawing figures..');
tic
if grouped_YN
    count_gene= histRead_group(allbt,exp_tags,exp_groups,...
        expected_list,'before QT',homomer_index_allbt);
    count_group = histRead(allbt,exp_groups,...
        expected_list,'before QT',homomer_index_allbt,1);
else
    count_gene = histRead(allbt,exp_tags,...
        expected_list,'before QT',homomer_index_allbt,1);
end
[thres,count] = barQuality_2(0,1,catagory_read_allbt,seq_quality_min_allbt,1);
drawnow;
toc

writefiles;

%% functions
% function to check input arguments, and import data
    function checkinputandimport
        if exist(input_file,'file')
        else error('Could not find the input file.');
        end
        
        switch length(varargin)
            case 0
                General_Alignment_ACGT_column_number = [4,5,6,7,8,9];
                XYPosition_ParentCell_column_number = [10,11,12];
            case 1
                error('Something is wrong. Contact Xiaoyan');
            case 2
                General_Alignment_ACGT_column_number = varargin{1};
                XYPosition_ParentCell_column_number = varargin{2};
        end
        
        if calculate_global_position_YN
            if (retiled_by_Python_or_MATLAB_YN && retiled_by_ZEN_YN)
                error('Specify how the images are tiled. Conflicting input arguments are detected.');
            elseif (retiled_by_Python_or_MATLAB_YN==0 && retiled_by_ZEN_YN==0)
                error('Specify how the images are tiled. Not enough input arguments.');
            else 
                if exist(csv_file_contain_tile_position,'file')
                else error('Could not find the position file.');
                end
            end
        end

        if cycle5_empty_YN
            if num_hybs==5
            else 
                error(['Uncompatible number of hybridization steps. '...
                    'If the 5th empty cycle exists, num_hybs has to be 5.']);
            end
        end
        
        % import data from files
        %--------------------
        seqdata = decodinginput(input_file,General_Alignment_ACGT_column_number,...
            XYPosition_ParentCell_column_number);
        
        if (calculate_global_position_YN && retiled_by_Python_or_MATLAB_YN)
            tile_start_pos = getposition(csv_file_contain_tile_position);
        end
        
        num_blobs = size(seqdata,1)/num_hybs;
        tile_ID_all = seqdata(:,3); 
        num_tiles = max(tile_ID_all);
        if floor(num_blobs)==num_blobs
        else
            error('Oops! The number of blobs is not integer. Check num_hybs again.');
        end
    end

% function to get quality items
    function qualityitems
        seq_quality(isnan(seq_quality)) = 0;

        if abnormal_sequencing_YN
            [seq_quality_min,seq_quality_min_ind] = min(seq_quality(:,inuni),[],2);
            [general_strength_min,general_strength_min_ind] = min(general_strength(:,inuni),[],2);
            [alignment_score_min,alignment_score_min_ind] = min(alignment_score(:,inuni),[],2);    
            if cycle5_empty_YN
                seq_res_b5 = seq_res(:,5);
                seq_quality_b5 = seq_quality(:,5);
                general_strength_b5 = general_strength(:,5);
                channel_max_b5 = channel_strength_max(:,5);
            end
        else
            if cycle5_empty_YN
                [seq_quality_min,seq_quality_min_ind] = min(seq_quality(:,1:4),[],2);
                [general_strength_min,general_strength_min_ind] = min(general_strength(:,1:4),[],2);
                [alignment_score_min,alignment_score_min_ind] = min(alignment_score(:,1:4),[],2); 
                seq_res_b5 = seq_res(:,5);
                seq_quality_b5 = seq_quality(:,5);
                general_strength_b5 = general_strength(:,5);
                channel_max_b5 = channel_strength_max(:,5);
            else
                [seq_quality_min,seq_quality_min_ind] = min(seq_quality,[],2);
                [general_strength_min,general_strength_min_ind] = min(general_strength,[],2);
                [alignment_score_min,alignment_score_min_ind] = min(alignment_score,[],2);
            end
        end
    end

% function to remove nonsensical reads
    function qualitycheck
        if check_parent_cell_YN && check_alignment_YN
            ffbt = (seq_quality_min>0 & general_strength_min>0 & cell_ID>0 & alignment_score_min>=alignment_min_threshold);
        elseif check_parent_cell_YN && check_alignment_YN==0
            ffbt = (seq_quality_min>0 & general_strength_min>0 & cell_ID>0);
        elseif check_parent_cell_YN==0 && check_alignment_YN
            ffbt = (seq_quality_min>0 & general_strength_min>0 & alignment_score_min>=alignment_min_threshold);
        else 
            ffbt = (seq_quality_min>0 & general_strength_min>0);
        end

        if nnz(ffbt)
        else
            error('No reads available. Try without parent cell check or alignment score.');
        end
    end

% function to write output files
    function writefiles
        disp('writing files..');
        tic
        output_directory_decode = [output_directory_decode '\'];
        if exist(output_directory_decode, 'dir')  
        else mkdir (output_directory_decode);
        end 

        % write code_n_count file
        %--------------------
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_code_n_count' '.csv'], 'w');
        if grouped_YN
            fprintf(fid,'Code,Count,GeneName,Group\n');
            temp_write = [cellstr(uni_letters_allbt),...
                num2cell(count_code_allbt'),...
                uni_reads_name_tag_allbt,...
                uni_reads_group_tag_allbt];
             for row = 1:size(temp_write,1)
                fprintf(fid,'%s,%d,%s,%s\n',temp_write{row,:});
            end
            
        else
            fprintf(fid,'Code,Count,GeneName\n');
            temp_write = [cellstr(uni_letters_allbt),...
                num2cell(count_code_allbt'),...
                uni_reads_name_tag_allbt];
            for row = 1:size(temp_write,1)
                fprintf(fid,'%s,%d,%s\n',temp_write{row,:});
            end
        end
        fclose(fid);

        % write gene_n_count file
        %--------------------
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_gene_n_count' '.csv'], 'w');
        if grouped_YN
            fprintf(fid,'GeneName,Count,Group\n');
            for row=1:length(count_gene(:,1))
                fprintf(fid,'%s,%d,%s\n',count_gene{row,:});
            end
        else
            fprintf(fid,'GeneName,Count\n');
            for row=1:length(count_gene(:,1))
                fprintf(fid,'%s,%d\n',count_gene{row,:});
            end
        end
        fclose(fid);

        % write group_n_count file
        %--------------------
        if grouped_YN
            fid = fopen([output_directory_decode output_filename_decode_prefix...
                '_group_n_count' '.csv'], 'w');
            fprintf(fid,'GeneGroup,Count\n');
            for row=1:length(count_group(:,1))
                fprintf(fid,'%s,%d\n',count_group{row,:});
            end
            fclose(fid);
        end

        % write qualitybar file
        %--------------------
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_qualitybar' '.csv'], 'w');
        fprintf(fid,'threshold,expected,unexpected,homomer,belowQT\n');
        temp_write = [thres',count];
        fprintf(fid,'%.2f,%d,%d,%d,%d\n',temp_write');
        fclose(fid);
        
        % write guessT file
        %--------------------
        if remove_T_in_analysis
            fid = fopen([output_directory_decode output_filename_decode_prefix...
                '_guessT' '.csv'], 'w');
            fprintf(fid,'%s\n','This is a table showing the guessing of T');
            fprintf(fid,'');
            for j = 1:length(Tcolumnn)
                fprintf(fid,',%s',Tcolumnn{j});
            end
            fprintf(fid,'\n');
            for i = 1:length(Trown)
                fprintf(fid,'%s,%d',Trown{i},Tdata{i,1});
                 for k = 2:length(Tcolumnn)
                    fprintf(fid,',%s',Tdata{i,k});
                end
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
        % write details file
        %--------------------
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_details' '.csv'], 'w');
        if calculate_global_position_YN
            if check_alignment_YN
                fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
                    'letters','name','global_X_pos','global_Y_pos',...
                    'parent_cell','tile_ID','general_stain_min',...
                    'seq_quality_min','alignment_score_min');
                temp_write = [cellstr(letters_allbt),name_tag_allbt,...
                    num2cell([global_x_allbt,global_y_allbt,...
                    cell_allbt,tile_allbt...
                    general_strength_min_allbt,seq_quality_min_allbt,...
                    alignment_score_min_allbt])];
                for row = 1:size(temp_write,1)
                    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d,%d\n',temp_write{row,:});
                end
                
            else
                fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n',...
                    'letters','name','global_X_pos','global_Y_pos',...
                    'parent_cell','tile_ID','general_stain_min',...
                    'seq_quality_min');
                temp_write = [cellstr(letters_allbt),name_tag_allbt,...
                    num2cell([global_x_allbt,global_y_allbt,...
                    cell_allbt,tile_allbt...
                    general_strength_min_allbt,seq_quality_min_allbt])];
                for row = 1:size(temp_write,1)
                    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d\n',temp_write{row,:});
                end
            end
        else
            if check_alignment_YN
                fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
                    'letters','name','tile_X_pos','tile_Y_pos',...
                    'parent_cell','tile_ID','general_stain_min',...
                    'seq_quality_min','alignment_score_min');
                temp_write = [cellstr(letters_allbt),name_tag_allbt,...
                    num2cell([tile_x_allbt,tile_y_allbt,...
                    cell_allbt,tile_allbt...
                    general_strength_min_allbt,seq_quality_min_allbt,...
                    alignment_score_min_allbt])];
                for row = 1:size(temp_write,1)
                    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d,%d\n',temp_write{row,:});
                end
            else
                fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n',...
                    'letters','name','tile_X_pos','tile_Y_pos',...
                    'parent_cell','tile_ID','general_stain_min',...
                    'seq_quality_min');
                temp_write = [cellstr(letters_allbt),name_tag_allbt,...
                    num2cell([tile_x_allbt,tile_y_allbt,...
                    cell_allbt,tile_allbt...
                    general_strength_min_allbt,seq_quality_min_allbt])];
                
                for row = 1:size(temp_write,1)
                    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d\n',temp_write{row,:});
                end
            end
            fclose(fid);
        end
        toc
        clear -regexp ^temp 
    end

%% save mat file
total_t2 = clock;
total_t = etime(total_t2,total_t1);
fprintf('%s%6f%s\n','Total elapsed time: ',total_t,'s');
clear -regexp ^temp ^total
disp('saving workspace variables..');
save([output_directory_decode output_filename_decode_prefix '.mat']); 

fprintf('Decoding beforeQT finished.\n\n');
end