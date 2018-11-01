function DecodingFunctionLite_Nbases(CpResultFile,TilePosFile,taglist,num_hybs,num_bases,quality_threshold)

% Light version of in situ sequencing decoding (function)
% More than four different types of bases are allowed.
% If one barcod is shared by multiple genes, only the first gene in the
% taglist will be assigend to all corresponding reads.
% No general stain thresholding.
% Xiaoyan, 2015-12-7


%% load files and extract info
disp('reading files');
seqdata = csvread(CpResultFile,1);
fID = fopen(TilePosFile);
tile_pos_file = textscan(fID, '%s%s%s%*[^\n]','headerlines',1,'delimiter',',');
tile_id = tile_pos_file{1};
xstart = tile_pos_file{2};
ystart = tile_pos_file{3};
for l = 1:length(tile_id)
    tile_start_pos(l,1) = str2double(tile_id{l});
    tile_start_pos(l,2) = str2double(xstart{l});
    tile_start_pos(l,3) = str2double(ystart{l});
end
fclose(fID);

%% extract info from taglist
for i = 1:size(taglist,1)
    exp_tags(i) = taglist(i,1);
    [token,remain] = strtok(exp_tags(i));
    exp_letters(i) = token;
    token = strtok(remain);
    exp_tags(i) = token;
end
expected_digits = letter2num(exp_letters,num_hybs); % child function letter2num in the end of the file
expected_list = expected_digits(:,1);

%% preallocate
num_blobs = size(seqdata,1)/num_hybs;
num_tiles = max(seqdata(:,3));

seq_res = zeros(num_blobs,num_hybs); 
channel_strength_max = zeros(num_blobs,num_hybs);
general_strength = zeros(num_blobs,num_hybs);
channel_strength_sum = zeros(num_blobs,num_hybs);
alignment_score = zeros(num_blobs,num_hybs);
tile_ID = zeros(num_blobs,1);
cell_ID = zeros(num_blobs,1);
global_x_pos = zeros(num_blobs,1);
global_y_pos = zeros(num_blobs,1);

%% DECODING
start_blob_ID=0;
fprintf('number of tiles\t number of blobs\n'); 
fprintf('\t%6d\t\t%6d\n',num_tiles,num_blobs);
fprintf('# tile\tnumber of blobs\n');

for t=1:num_tiles
    temp_tile_seq_data = seqdata(seqdata(:,3) == t,:);  % blobs within a tile

    temp_tile_num_blobs = length(temp_tile_seq_data(:,2))/num_hybs;

    temp_tile_x_start = mean(tile_start_pos(tile_start_pos(:,1) == t,2));
    temp_tile_y_start = mean(tile_start_pos(tile_start_pos(:,1) == t,3));
   
    fprintf('%6d\t\t%6d\n',t,temp_tile_num_blobs);

    %------------------------
    if temp_tile_num_blobs
        [~,order] = sort(temp_tile_seq_data(:,2));  % order based on the object number
       
        % extract info about blobs within the tile
        temp_tile_intensity = temp_tile_seq_data(order,6:5+num_bases);
        temp_tile_general = temp_tile_seq_data(order,4);
        temp_tile_align = temp_tile_seq_data(order,5);
        temp_tile_x_pos = temp_tile_seq_data(order,6+num_bases);
        temp_tile_y_pos = temp_tile_seq_data(order,7+num_bases);
        
                
        % lower T intensity of base3
        temp_tile_intensity((temp_tile_num_blobs*2+1):temp_tile_num_blobs*3, 4) = ...
            temp_tile_intensity((temp_tile_num_blobs*2+1):temp_tile_num_blobs*3, 4)/3;
        
        % max intensities
        [maxint,maxidx] = max(temp_tile_intensity,[],2);
        channel_strength_max(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(maxint,num_hybs,temp_tile_num_blobs))';

        % seq result matrix
        seq_res(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(maxidx,num_hybs,temp_tile_num_blobs))';
        
        % blob info
        tile_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = t;
        alignment_score(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(temp_tile_align,num_hybs,temp_tile_num_blobs))';
        
        temp_list = 1:temp_tile_num_blobs*num_hybs;
        if num_hybs ~= 1
            temp_list = mod(temp_list,num_hybs)==1;  % take the information only from the first base
        end
        
        global_x_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            repmat(temp_tile_x_start,temp_tile_num_blobs,1) + temp_tile_x_pos(temp_list);
        global_y_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            repmat(temp_tile_y_start,temp_tile_num_blobs,1) + temp_tile_y_pos(temp_list);
        
        general_strength(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
           (reshape(temp_tile_general,num_hybs,temp_tile_num_blobs))';
        channel_strength_sum(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            reshape(sum(temp_tile_intensity,2),num_hybs,temp_tile_num_blobs)';
        
        start_blob_ID = start_blob_ID+temp_tile_num_blobs;
    end

end
clear -regexp ^temp
clear t c

%% combine into reads
seq_num = zeros(num_blobs,1);
j = num_hybs;
for i = 1:num_hybs
    seq_num = seq_res(:,i).*10^(j-1)+seq_num;
    j = j-1;
end

%% sequencing quality
seq_quality = channel_strength_max./channel_strength_sum;
seq_quality = min(seq_quality,[],2);
general_strength = min(general_strength,[],2);

%% remove nonsensical reads
ffbt = (seq_quality>0 & general_strength>0);    
allbt = seq_num(ffbt); 
seq_res_allbt = seq_res(ffbt,:);
letters_allbt = num2letters(seq_res_allbt); % child function num2letters in the end of the file
general_strength_min_allbt = general_strength(ffbt,:);
seq_quality_min_allbt = seq_quality(ffbt,:);
tile_allbt = tile_ID(ffbt);
global_x_allbt = global_x_pos(ffbt);
global_y_allbt = global_y_pos(ffbt);

%% assign gene names
name_tag_allbt = cell(length(allbt),1);
name_tag_allbt(~ismember(allbt,expected_list)) = {'NNNN'};
for i = 1:length(expected_list)
    name_tag_allbt(allbt == expected_list(i)) = exp_tags(i);
end

%% figure - before QT
[~,count_code_allbt,uni_letters_allbt,uni_reads_name_tag_allbt,~,catagory_read_allbt] = ...
    groupreads(allbt,seq_res_allbt,letters_allbt,expected_list,exp_tags,cell(length(exp_tags),1));

homomer_index_allbt = findHomomer(expected_list,allbt);
catagory_read_allbt(homomer_index_allbt) = -1;

count_gene = histRead(allbt,exp_tags,expected_list,'before QT',homomer_index_allbt,1);
[~,~] = barQuality_2(0,1,catagory_read_allbt,seq_quality_min_allbt,1);
drawnow;

%% write files - before QT
disp('writing files - before QT')

fid = fopen('code_n_count_beforeQT.csv','w');
fprintf(fid,'Code,Count,GeneName\n');
temp_write = [cellstr(uni_letters_allbt),...
    num2cell(count_code_allbt'),...
    uni_reads_name_tag_allbt];
for row = 1:size(temp_write,1)
    fprintf(fid,'%s,%d,%s\n',temp_write{row,:});
end
fclose(fid);

fid = fopen('gene_n_count_beforeQT.csv','w');
fprintf(fid,'GeneName,Count\n');
for row=1:length(count_gene(:,1))
    fprintf(fid,'%s,%d\n',count_gene{row,:});
end
fclose(fid);

fid = fopen('details_beforeQT.csv','w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',...
    'letters','name','global_X_pos','global_Y_pos',...
    'tile_ID','general_stain_min',...
    'seq_quality_min');
temp_write = [cellstr(letters_allbt),name_tag_allbt,...
    num2cell([global_x_allbt,global_y_allbt,...
    tile_allbt...
    general_strength_min_allbt,seq_quality_min_allbt])];
for row = 1:size(temp_write,1)
    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d\n',temp_write{row,:});
end
fclose(fid);

%% quality thresholding
ff = (seq_quality_min_allbt>quality_threshold);

if ~nnz(ff)
    error('No reads available. Try a lower threshold.');
end

allqt = allbt(ff);
seq_res_allqt = seq_res_allbt(ff,:);
letters_allqt = letters_allbt(ff,:);
name_tag_allqt = name_tag_allbt(ff);

seq_quality_min_allqt = seq_quality_min_allbt(ff,:);
general_strength_min_allqt = general_strength_min_allbt(ff,:);

tile_allqt = tile_allbt(ff);
global_x_allqt = global_x_allbt(ff);
global_y_allqt = global_y_allbt(ff);

%% figure - after QT
[~,count_code_allqt,uni_letters_allqt,uni_reads_name_tag_allqt,~,~] = ...
    groupreads(allqt,seq_res_allqt,letters_allqt,expected_list,exp_tags,cell(length(exp_tags),1));

homomer_index_allqt = findHomomer(expected_list,allqt);
count_gene = histRead(allqt,exp_tags,expected_list,['QT - ' num2str(quality_threshold)],homomer_index_allqt,1);

%% write files - after QT
disp('writing files - after QT')

fid = fopen(['code_n_count_QT_' num2str(quality_threshold) '.csv'],'w');
fprintf(fid,'Code,Count,GeneName\n');
temp_write = [cellstr(uni_letters_allqt),...
    num2cell(count_code_allqt'),...
    uni_reads_name_tag_allqt];
for row = 1:size(temp_write,1)
    fprintf(fid,'%s,%d,%s\n',temp_write{row,:});
end
fclose(fid);

fid = fopen(['gene_n_count_QT_' num2str(quality_threshold) '.csv'],'w');
fprintf(fid,'GeneName,Count\n');
for row=1:length(count_gene(:,1))
    fprintf(fid,'%s,%d\n',count_gene{row,:});
end
fclose(fid);

fid = fopen(['details_QT_' num2str(quality_threshold) '.csv'],'w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',...
    'letters','name','global_X_pos','global_Y_pos',...
    'tile_ID','general_stain_min',...
    'seq_quality_min');
temp_write = [cellstr(letters_allqt),name_tag_allqt,...
    num2cell([global_x_allqt,global_y_allqt,...
    tile_allqt...
    general_strength_min_allqt,seq_quality_min_allqt])];
for row = 1:size(temp_write,1)
    fprintf(fid,'%s,%s,%d,%d,%d,%d,%d\n',temp_write{row,:});
end
fclose(fid);

%% child functions
function out=letter2num(letter_list,num_hybs)

letters = {'A';'C';'G';'T';'K'};
tag_list = [];
for i = 1:length(letter_list)
    tag = []; 
    tag_temp = letter_list(i);
    tag_temp = tag_temp{1};
    
     for j = 1:num_hybs
        tag_lett = tag_temp(j);
        a = find(strcmp(letters,tag_lett));
        
        if a
            taglist_num(i,j) = a;
            a = num2str(a);
            tag = [tag a];
        else 
            error('Unexpected barcode letter found.');
        end
       
     end
     
     tag_list = [tag_list; tag]; 
end

exp_list = str2num(tag_list);

out = [exp_list taglist_num];

end


function m_letters = num2letters(m_num)

letters = ['A';'C';'G';'T';'K';'O'];
m_letters = repmat(char(0),size(m_num,1),size(m_num,2));

for i = 1:size(m_num,1)
    for b=1:size(m_num,2)
        m_letters(i,b) = letters(m_num(i,b));
    end
end
end

end