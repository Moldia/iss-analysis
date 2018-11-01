fixdata1 = importdata('CP_160919_fix\blobs.csv');
fixdata1 = fixdata1.data;

fixdata2 = importdata('CP_160920_fix2\blobs.csv');
fixdata2 = fixdata2.data;

data = importdata('blobs.csv');
headers = data.textdata;
data = data.data;

%%
for i = 45:62
    % base1
    data(data(:,1)==(i-1)*4+1,:) = fixdata1(fixdata1(:,3)==i & mod(fixdata1(:,1),2)==1,:);
    % base3
    data(data(:,1)==(i-1)*4+3,:) = fixdata1(fixdata1(:,3)==i & mod(fixdata1(:,1),2)==0,:);
end

for i = 128:162
    data(data(:,1)==(i-1)*4+1,:) = fixdata2(fixdata2(:,3)==i,:);
end

for i = 169:201
    data(data(:,1)==(i-1)*4+3,:) = fixdata2(fixdata2(:,3)==i,:);
end



%%
headers = strcat(headers,',');
headers = [headers{:}];
headers = [headers(1:end-1), '\n'];

fid = fopen('blobs_fixed.csv', 'w');
fprintf(fid,headers);
fmt = repmat('%d,',1,size(data,2));
fmt = [fmt(1:end-1), '\n'];
fprintf(fid,fmt,data');
fclose(fid);