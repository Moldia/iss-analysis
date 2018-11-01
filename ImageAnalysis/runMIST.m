ntilesX = 17;
ntilesY = 12;
emptytiles = csvread('EmptyTiles.csv', 1);
emptytiles = ntilesX*(emptytiles(:,2)-1)+emptytiles(:,1);
imwrite(zeros(2048, 'uint16'), 'FS_tophat_stack\empty.tif');

img_name_grid = cell(ntilesX*ntilesY, 1);
for i = 1:ntilesX*ntilesY
    if ismember(i, emptytiles)
        img_name_grid{i} = 'empty.tif';
    else
        % tile images
        img_name_grid{i} = ['Ab_t', paddigits(i-nnz(emptytiles<i),3), '_FS_stack.tif'];
    end
end
img_name_grid = reshape(img_name_grid, ntilesX, ntilesY)';


% % use MIST's original function to build image grid
% % all files in the grid need to exist
% [img_name_grid, to_stitch_time_slice_nbs] = build_img_name_grid...
%     ('K:\161101_4028_36-2\MIP_tophat', 'base1_c1t{ppp}_MIP.tif',...
%     ntilesX, ntielsY,...
%     'tiling_technique','combing',...
%     'starting_point','upperleft',...
%     'first_direction','horizontal');

        
input_directory = 'FS_tophat_stack';
output_directory = 'Stitched';
output_prefix = 'stitchtest';
save_stitched_image = 1;
assemble_from_metadata = 0;     % apply the same to other images

% % no overlap
% create_zero_percent_overlap_mosaic(input_directory, img_name_grid,...
%     output_directory, output_prefix, 'Max', 0, 1);

stitch_time_slice(input_directory, img_name_grid,...
    output_directory, output_prefix,...
    1, NaN, NaN,...     % these will probably not change
    'Max', 0,...        % blending method
    save_stitched_image, assemble_from_metadata,...
    [output_directory '\logtest.txt'],...
    10, 10)

