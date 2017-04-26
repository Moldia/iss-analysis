% tranform reference image and spot coordinates based on given rotation and
% transformation parameters
% Xiaoyan, 2017


clear;
close all;

%% input
image_to_transform = '..\Preprocessing\Stitched\Ch026_SBL2_c1_stitched.tif';
coordinates_file_to_transform = 'Ch026\Decoding\QT_0.5_0.0001_details.csv';

scaling = 1;
rotation = 180;     % in degree, positive: counterclockwise
translation = [0, 0];     % shift in x and y, positive: left (for x) and up (for y)

output_folder = 'Transformed';


%%
% organize data
I = imread(image_to_transform);
I = imresize(I, scaling);
imsize = size(I);

% transform image
mkdir(output_folder);
name = strsplit(image_to_transform, filesep);
name = strsplit(name{end}, '.tif');
name = fullfile(output_folder, [name{1}, '_Transformed.tif']);
if numel(size(I)) > 2
    I2 = [];
    for j = 1:3
        Itemp = transformimage(I(:,:,j),...
            rotation, translation(2), translation(1), 1, name, imsize, 0);
        I2(:,:,j) = Itemp;
    end
    imwrite(I2, name);
else
    I2 = transformimage(I, rotation, translation(2), translation(1), 1, name);
end

% transform coordinates
[~, pos] = getinsitudata(coordinates_file_to_transform);
pos2 = transform_coordinates(pos, scaling, imsize, rotation,...
    translation(2), translation(1), size(I2), 1);

% save file
name = strsplit(coordinates_file_to_transform, filesep);
name = strsplit(name{end}, '.csv');
name = fullfile(output_folder, [name{1}, '_Transformed.csv']);
replace_file_columns(coordinates_file_to_transform, [3,4], [{pos(:,1)},{pos(:,2)}], name);

% visualization
figure; 
imshow(I2, []); hold on;
plot(pos2(:,1), pos2(:,2), '.');


