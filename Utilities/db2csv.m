clear

load 417

Name = {'Aligned_SpecBlob' 'AlignmentMath' 'aligned_A' 'aligned_C' 'aligned_G'};

temp_intensity = {};
for j = 1:5
    field = ['Intensity_MaxIntensity_' Name{j}];
    temp_intensity = [temp_intensity,(getfield(handles.Measurements.blobs,field))'];
end
% temp_intensity = cell2mat(temp_intensity);

temp_location = [(handles.Measurements.blobs.Location_Center_X)',...
    (handles.Measurements.blobs.Location_Center_Y)'];
temp_parent = handles.Measurements.blobs.Parent_Cells;
temp_meata = handles.Measurements.Image.Metadata_position;

Write = double([]);
for i = 1:416
    temp_write = double(temp_parent{i});
    Write = [Write;repmat(i,size(temp_write,1),1),(1:size(temp_write,1))',...
        repmat(double(temp_meata{i}),size(temp_write,1),1),...
        cell2mat(temp_intensity(i,:)),...
        cell2mat(temp_location(i,:)),...
        temp_write];
end

Title = [{'ImageNumber' 'ObjectNumber' 'Metadata_position'} Name {'Location_Center_X' 'Location_Center_Y' 'Parent_cell'}];

fid = fopen('blobs_merged.csv','w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',Title{:});
fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',Write');
fclose(fid);
