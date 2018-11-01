% combine multiple CP position input files
% all of them must have the same header and column numbers
% Xiaoyan, 2017


%% input
files = {'E:\PROOOJECTS\1_Neuron_mapping\Experiments\160408_CorticalProbes\160222UCL_GCaMP-hippo3\ZENout\Tiled_160417.csv',...
    'E:\PROOOJECTS\1_Neuron_mapping\Experiments\160408_CorticalProbes\160222UCL_GCaMP-hippo3\ZENout\Tiled_160417 - Copy.csv'};
startingpos_x = [0, 1000, 2000];
startingpos_y = [0, 20000, 40000];

%%
Data = [];
iTileStart = 0;

for i = 1:length(files)
    % read files and keep original header
    data = importdata(files{i});  % normally will all be recognized as text
    data = cellfun(@(v) strsplit(v, ','), data, 'uni', 0);
    data = {cat(1, data{:})};
    header = data{1}(1,:);
    data = data{1}(2:end,:);    % data: column1-tile#, column2-xpos, column3-ypos
    
    % remove extra tiles which resulted from padding or incorrect exporting 
    data(...
        startingpos_x(i)+cellfun(@str2double, data(:,2)) >= startingpos_x(i+1)...
        & ...
        startingpos_y(i)+cellfun(@str2double, data(:,3)) >= startingpos_y(i+1),:) = [];

    % tile #
    data(:,1) = num2cell(iTileStart + cellfun(@str2double, data(:,1)));
    % start pos x
    data(:,2) = num2cell(startingpos_x(i) + cellfun(@str2double, data(:,2)));
    % start pos y
    data(:,3) = num2cell(startingpos_y(i) + cellfun(@str2double, data(:,3)));
    
    Data = [Data; data];
    iTileStart = max(cell2mat(Data(:,1)));
end

%% write file
Data = Data';
fmt = '%d,%d,%d';
for i = 1:length(header)-3
    fmt = [fmt ',%s'];
end
header = strcat(header, ',');
header = [header{:}];

fid = fopen('Tiled_comb.csv', 'w');
fprintf(fid, [header(1:end-1) '\n']);
fprintf(fid, [fmt '\n'], Data{:});
fclose(fid);
