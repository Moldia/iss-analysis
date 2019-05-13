%% files to read
samples = ls('..');
samples =  cellstr(samples(3:14,:));

%% read and import
All = cell(length(samples),3);
for i = 1:length(samples)
    decodedir = strcat('..\',samples{i},'\Decoding\');
    decodefile = strcat(decodedir, 'QT_0.4_0.005_details.csv');
    [name,pos,quality] = getinsitudata_f(decodefile);
    All(i,:) = [{name},{pos},{quality}];
end

%% remove NNNN, bin data and count transcripts
All_re = cell(length(samples),3);
Name_all = {};
Name_map = cell(length(samples),1);
grid_size = 600;
Count = cell(length(samples),1);
% figure;
for i = 1:length(samples)
    i
    [name_uni,~,idx_name] = unique(All{i,1});
%     subplot(3,4,i); hold on; axis image;
%     for j = 1:length(name_uni)
%         plot(All{i,2}(idx_name==j,1),All{i,2}(idx_name==j,2),'+');
%     end
    idx_NNNN = find(strcmp(name_uni,'NNNN'));
%     subplot(3,4,i);
%     plot(All{i,2}(idx_name==idx_NNNN,1),All{i,2}(idx_name==idx_NNNN,2),'+')
    
    noNNNN = cellfun(@(v) v(idx_name~=idx_NNNN,:),All(i,:),'uni',0);    % remove NNNN
    All_re(i,:) = noNNNN;
    name = noNNNN{1};
    pos = noNNNN{2};
    quality = noNNNN{3};
    
    [name_uni,~,idx_name] = unique(name);
    
    % collect all transcripts and map between indeces
    Name_all = [Name_all; name_uni(~ismember(name_uni,Name_all))];
    temp = cellfun(@(v) strcmp(v,Name_all),name_uni,'uni',0);
    Name_map{i} = cellfun(@find,temp);
    
    % divide into grids
    grid_nrx =  ceil(max(pos(:,1))/grid_size);
    grid_nry =  ceil(max(pos(:,2))/grid_size);
    grid_nr = grid_nrx*grid_nry;

    Pos_pixel = round(correctcoordinates_f(pos,1/grid_size));

    blobinpolygon = cell(1,grid_nr);
    Counted = false(length(pos),1);    
    for j = 1:grid_nrx   % j along x axis, column
%         fprintf('%.2f%s\n',double((j-1)/grid_nrx)*100,'% grid processed.');
        for k = 1:grid_nry  % i along y axis, row
            temp_in = Pos_pixel(:,1)==j & Pos_pixel(:,2)==k;
            temp_in = temp_in & ~Counted;
            
            % blobs within grid (logical)
            p = (j-1)*grid_nry+k; % counting direction: y
            blobinpolygon(p) = {temp_in};
            
            % already counted blobs
            Counted = Counted | temp_in;
            
        end
    end
    disp('100.00% grid processed.');
    
    % count transcripts in the grid
    count_transcript = zeros(length(name_uni),grid_nr);
    % tic
    for j = 1:grid_nr
        count_transcript(:,j) = (hist(idx_name(logical(blobinpolygon{j})),1:length(name_uni)))';
    end
    % toc
    
    emptygrid = sum(count_transcript,1)==0;
    count_transcript_nozero = count_transcript(:,~emptygrid);
    Count{i} = count_transcript_nozero;

end

%% remap transcripts
allgrids = cellfun(@(v) size(v,2),Count);
Count_all = zeros(sum(allgrids),length(Name_all)+1);

row = 0;
for i = 1:length(samples)
    nmap = Name_map{i};
    Count_all(row+1:row+allgrids(i),nmap) = Count{i}';
    Count_all(row+1:row+allgrids(i),end) = i;
    row = row + allgrids(i);
end
hist(sum(Count_all(:,1:end-1),2))

Count_norm = Count_all(sum(Count_all(:,1:end-1),2)>=10,:);
Count_norm(:,1:end-1) = bsxfun(@rdivide,Count_norm(:,1:end-1),sum(Count_norm(:,1:end-1),2));

%% reorganize and write files
headers = {'name','pos_x','pos_y','quality'};
for i = 1:length(samples)
    nmap = Name_map{i};
    names = All_re{i,1};
    [name_uni,~,idx_name] = unique(names);
    names_new = num2cell(idx_name);
    names_new = cellfun(@(v) ['gene' num2str(nmap(v))],names_new,'uni',0);
    towrite = [names_new,num2cell(All_re{i,2}),num2cell(All_re{i,3})]';
    fid = fopen([samples{i}, '_insitu_details.csv'],'w');
    fprintf(fid,'%s,%s,%s,%s\n',headers{:});
    fprintf(fid,'%s,%d,%d,%d\n',towrite{:});
    fclose(fid);
end

fid = fopen('All_genes.csv','w');
fprintf(fid,'%s\n',Name_all{:});
fclose(fid);

%% tSNE
mappedtsne = tsne(Count_norm(:,1:end-1), Count_norm(:,end), 2);
figure,gscatter(mappedtsne(:,1),mappedtsne(:,2),Count_norm(:,end));
legend(samples);
save('tsne.mat','Count_all','Count_norm','mappedtsne','Name_all');

mappedtsne_sub = tsne(Count_norm(Count_norm(:,end)==9|Count_norm(:,end)==10,1:end-1), Count_norm(Count_norm(:,end)==9|Count_norm(:,end)==10,end), 2);
