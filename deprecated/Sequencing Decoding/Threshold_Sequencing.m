function Threshold_Sequencing(output_directory_decode,...
    output_filename_decode_prefix,output_filename_afterQT_prefix,...
    quality_threshold,general_strain_threshold,remove_homomer_YN,...
    based_on_empty_cycle_5_YN,remove_only_above_quality_threshold_YN,varargin)

% Take reads above threshold.
% Two thresholds: general stain inteisty and sequencing quality.
% Can remove homomer reads.
%
% Option: remove homomer reads
%   SubOption: based on empty cycle 5.
%       Remove only the reads show the same base also in the 5th empty
%       cycle. The original experiment should contain this information and
%       also in decoding this should have been specified.
%   SubSubOption: remove only above quality threshold.
%       In order for a read to be removed in this way, it should show as 
%       5-base homomer. And in the empty cycle, it has a quality score
%       higher than the set threshold.
%
% Thershold_Sequencing v2.4
% Xiaoyan, 2014-11-29

%% initiate
drawnow;
disp('Initiating Threshold_Sequencing.');
output_directory_decode = [output_directory_decode '\'];
if exist([output_directory_decode output_filename_decode_prefix '.mat'],'file')~=2
    error(['Could not find the beforeQT file. '...
        'Make sure the working directory is correct.']);
else
    load([output_directory_decode output_filename_decode_prefix '.mat']);
end

%% quality threshold
ff = (seq_quality_min_allbt>quality_threshold & ...
    general_strength_min_allbt>general_strain_threshold);
count_with_homomer = nnz(ff);

% remove unexpected homomer reads
if remove_homomer_YN
    fprintf('Removing homomer reads.\n'); 
    if based_on_empty_cycle_5_YN
        if exist('seq_quality_b5','var')
            if remove_only_above_quality_threshold_YN
                f5 = (seq_quality_b5_bt(homomer_index_allbt)>quality_threshold...
                    & seq_res_allbt(homomer_index_allbt,1)==seq_res_b5_bt(homomer_index_allbt));
                ff(homomer_index_allbt(f5))=0;
            else
                f5 = (seq_res_allbt(homomer_index_allbt,1)==seq_res_b5_bt(homomer_index_allbt));
                ff(homomer_index_allbt(f5))=0;
            end
        else
            error(['No empty cycle5 information. '...
                'Check Decoding settings and run again from Decoding.']);
        end
    else
        ff(homomer_index_allbt)=0;
    end
    homomer_deleted = count_with_homomer - length(ff(ff==1));
    fprintf('count with homomer\t homomer deleted\n'); 
    fprintf('%10d \t\t\t %10d\n',[count_with_homomer homomer_deleted]);
end

if ~nnz(ff)
    error('No reads available. Try a lower threshold.');
end

%% all reads after QT
allqt = allbt(ff); 
seq_res_allqt = seq_res_allbt(ff,:);
letters_allqt = letters_allbt(ff,:);
name_tag_allqt = name_tag_allbt(ff);
group_tag_allqt = group_tag_allbt(ff);

seq_quality_allqt = seq_quality_allbt(ff,:);
general_strength_allqt = general_strength_allbt(ff,:);
alignment_score_allqt = alignment_score_allbt(ff,:);

channel_strength_max_allqt = channel_strength_max_allbt(ff,:);
channel_strength_original_allqt = channel_strength_original_allbt(:,:,ff);

tile_allqt = tile_allbt(ff);
tile_cell_allqt = tile_cell_allbt(ff);
cell_allqt = cell_allbt(ff);
if exist('tile_blob_allbt','var')
    tile_blob_allqt = tile_blob_allbt(ff);
end
if exist('blob_allbt','var')
    blob_allqt = blob_allbt(ff);
end
tile_x_allqt = tile_x_allbt(ff);
tile_y_allqt = tile_y_allbt(ff);
if calculate_global_position_YN
    global_x_allqt = global_x_allbt(ff);
    global_y_allqt = global_y_allbt(ff);
end

seq_quality_min_allqt = seq_quality_min_allbt(ff,:);
general_strength_min_allqt = general_strength_min_allbt(ff,:);
alignment_score_min_allqt = alignment_score_min_allbt(ff,:);

if exist('seq_quality_min_ind_allbt','var')
    seq_quality_min_ind_allqt = seq_quality_min_ind_allbt(ff,:);
end
if exist('general_strength_min_ind_allbt','var')
    general_strength_min_ind_allqt = general_strength_min_ind_allbt(ff,:);
end
if exist('alignment_score_min_ind_allbt','var')
    alignment_score_min_ind_allqt = alignment_score_min_ind_allbt(ff,:);
end

%% group reads
[uni_seq_res_allqt,count_code_allqt,uni_reads_letters_allqt,uni_reads_name_tag_allqt,uni_reads_group_tag_allqt,catagory_read_allqt] = ...
    groupreads(allqt,seq_res_allqt,letters_allqt,expected_list,exp_tags,exp_groups);

% find unexpected homomer reads
homomer_index_allqt = findHomomer(expected_list,allqt);
catagory_read_allqt(homomer_index_allqt) = -1;

%% output figures
disp('busy drawing figures..');
tic
if grouped_YN
    count_gene = histRead_group(allqt,exp_tags,exp_groups,expected_list,...
        ['QT-' num2str(quality_threshold)],homomer_index_allqt);
    count_group = histRead(allqt,exp_groups,expected_list,...
        ['QT-' num2str(quality_threshold)],homomer_index_allqt,1);
else
    count_gene = histRead(allqt,exp_tags,expected_list,...
        ['QT-' num2str(quality_threshold)],homomer_index_allqt,1);
end
drawnow;
toc

%% output files
disp('writing files..');
tic
% write code_n_count file
%--------------------
fid = fopen([output_directory_decode output_filename_afterQT_prefix...
    '_code_n_count' '.csv'], 'w');
if grouped_YN
    fprintf(fid,'Code,Count,GeneName,Group\n');
    temp_write = [cellstr(uni_reads_letters_allqt),...
                num2cell(count_code_allqt'),...
                uni_reads_name_tag_allqt,...
                uni_reads_group_tag_allqt];
    for row = 1:size(temp_write,1)
        fprintf(fid,'%s,%d,%s,%s\n',temp_write{row,:});
    end
else
    fprintf(fid,'Code,Count,GeneName\n');
    temp_write = [cellstr(uni_reads_letters_allqt),...
                num2cell(count_code_allqt'),...
                uni_reads_name_tag_allqt];
    for row = 1:size(temp_write,1)
        fprintf(fid,'%s,%d,%s\n',temp_write{row,:});
    end
end
fclose(fid);

% write gene_n_count file
%--------------------
fid = fopen([output_directory_decode output_filename_afterQT_prefix...
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
    fid = fopen([output_directory_decode...
        output_filename_afterQT_prefix '_group_n_count' '.csv'], 'w');
    fprintf(fid,'GeneGroup,Count\n');
    for row=1:length(count_group(:,1))
        fprintf(fid,'%s,%d\n',count_group{row,:});
    end
    fclose(fid);
end

% write details file
%--------------------
fid = fopen([output_directory_decode output_filename_afterQT_prefix...
    '_details' '.csv'], 'w');
if calculate_global_position_YN
    if check_alignment_YN
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
            'letters','name','global_X_pos','global_Y_pos',...
            'parent_cell','tile_ID','general_stain_min',...
            'seq_quality_min','alignment_score_min');
        temp_write = [cellstr(letters_allqt),name_tag_allqt,...
            num2cell([global_x_allqt,global_y_allqt,...
            cell_allqt,tile_allqt...
            general_strength_min_allqt,seq_quality_min_allqt,...
            alignment_score_min_allqt])];
        for row = 1:size(temp_write,1)
            fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d,%d\n',temp_write{row,:});
        end
    else
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n',...
            'letters','name','global_X_pos','global_Y_pos',...
            'parent_cell','tile_ID','general_stain_min',...
            'seq_quality_min');
        temp_write = [cellstr(letters_allqt),name_tag_allqt,...
            num2cell([global_x_allqt,global_y_allqt,...
            cell_allqt,tile_allqt...
            general_strength_min_allqt,seq_quality_min_allqt])];
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
        temp_write = [cellstr(letters_allqt),name_tag_allqt,...
            num2cell([tile_x_allqt,tile_y_allqt,...
            cell_allqt,tile_allqt...
            general_strength_min_allqt,seq_quality_min_allqt,...
            alignment_score_min_allqt])];
        for row = 1:size(temp_write,1)
            fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d,%d\n',temp_write{row,:});
        end
    else
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n',...
            'letters','name','tile_X_pos','tile_Y_pos',...
            'parent_cell','tile_ID','general_stain_min',...
            'seq_quality_min');
        temp_write = [cellstr(letters_allqt),name_tag_allqt,...
            num2cell([tile_x_allqt,tile_y_allqt,...
            cell_allqt,tile_allqt...
            general_strength_min_allqt,seq_quality_min_allqt])];
        for row = 1:size(temp_write,1)
            fprintf(fid,'%s,%s,%d,%d,%d,%d,%d,%d\n',temp_write{row,:});
        end
    end
end
fclose(fid);
toc

% save mat file
%--------------------
clear -regexp ^temp
clear fid i 
disp('saving workspace variables..');
save([output_directory_decode output_filename_afterQT_prefix '.mat']); 

fprintf('Thresholding finished.\n\n'); 
end

