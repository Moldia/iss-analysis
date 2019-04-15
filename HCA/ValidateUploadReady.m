codebook = 'G:\HCA_09b_human_120genes\starfish_codebook_human120.csv';
folder = 'G:\HCA_09b_human_120genes\upload_starfish_190415\ISS_h_brain_03\main_files';
tile = 345;  % example tile

%%
code = importdata(codebook, ',');
genes = code.textdata;
code = code.data;

% pixel calling
reads = [];
for r = 1:size(code,2)
    Ibase = [];
    for c = 1:max(code(:))
        I = imread(fullfile(folder, ['primary_images-fov_' paddigits(tile-1,3) '-Z0-H' num2str(r-1) '-C' num2str(c-1) '.tiff']));
        Ibase = cat(3, Ibase, imtophat(I, strel('disk', 3)));
    end
    [~, base] = max(Ibase, [], 3);
    reads = cat(3, reads, base);
end

Reads = permute(reads, [3 1 2]);
Reads = reshape(Reads, size(code,2), [])';

Genes = zeros(2000, 2000);
Expected = ismember(Reads, code, 'rows');
temp = find(Expected);
for i = 1:size(code,1)
    idx = ismember(Reads(temp,:), code(i,:), 'rows');
    [i, nnz(idx)]
    Genes(temp(idx)) = i;
    temp(idx) = [];
end
        

%% visualize most abundent code that has at least three colors
count = hist(Genes(:), 0:length(genes));
allfour = find(cellfun(@(v) length(unique(v)), num2cell(code, 2))>=3);
[~, idx] = sort(count(allfour+1), 'descend');
idx = allfour(idx(1));

figure; Ax = [];
[y,x] = ind2sub([2000, 2000], find(Genes == idx));
for r = 1:5
    ax = subplot(2,3,r); Ax = [Ax, ax];
    I = imread(fullfile(folder, ['primary_images-fov_' paddigits(tile-1,3) '-Z0-H' num2str(r-1) '-C' num2str(code(idx,r)-1) '.tiff']));
    imshow(I, []); hold on; plot(x,y,'.');
    title(['H' num2str(r-1) '-C' num2str(code(idx,r)-1)]);
end

% ax = subplot(236); Ax = [Ax, ax];
% imagesc(Genes); 
% colormap(gca, 'hot');
% ch = colorbar;
% ch.Ticks = 1:length(genes);
% ch.TickLabels = genes;
% axis image;
% axis tight;

ax = subplot(236); Ax = [Ax, ax];
I = imread(fullfile(folder, ['nuclei-fov_' paddigits(tile-1,3) '-Z0-H0-C0.tiff']));
imshow(I, []); hold on; plot(x,y,'.');
title(['auxiliary ' genes{idx}]);
    
linkaxes(Ax, 'xy');
        

