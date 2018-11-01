% combine multiple blob files which have the same header
% Xiaoyan, 2017

folders = {'CP_160917_1' 'CP_160917_2', 'CP_160917_3'};
files = {'CP_160917_1\blobs_fixed.csv' 'CP_160917_2\blobs.csv', 'CP_160917_3\blobs.csv'};

Data = [];
for i = 1:length(files)
%     data = importdata([folders{i}, '\blobs.csv']);
    data = importdata(files{i});
    headers = data.textdata;
    Data = [Data; data.data];
end


headers = strcat(headers, ',');
headers = [headers{:}];
headers = [headers(1:end-1), '\n'];

% write
fid = fopen('blobs_comb.csv', 'w');
fprintf(fid,headers);
fmt = repmat('%d,',1,size(Data,2));
fmt = [fmt(1:end-1), '\n'];
fprintf(fid,fmt,Data');
fclose(fid);
