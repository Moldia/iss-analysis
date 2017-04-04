% kernel density estimation of all genes
% Xiaoyan, 2017

clear;
close all;

%% parameters
decoded_file = '..\CombinedResults_wCell.csv';
image = '..\Preprocessing\Stitched\AlignedImages\base2_c1_ORG.tif';   % important for size
scale = .2;
bandwidth = 30; % in original scale


%%
% transcripts
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);

% unique transcripts
[uniName, ~, idxName] = unique(name);
countName = hist(idxName, 1:length(uniName));

% % remove reads with low counts
% [name, pos] = removereads(name, uniName(countName<200), pos);
% [uniName, ~, idxName] = unique(name);
% countName = hist(idxName, 1:length(uniName));

% image size
imsize = imfinfo(image);
imsize = [imsize.Height, imsize.Width];

% density estimation plot
pos = correctcoord(pos, .2);
f = 0;
figure;
set(gcf, 'units', 'normalized', 'position', [.05+.02*f .1-.02*f .8 .8]);
for i = 1:length(uniName)
    if mod(i-1,12)==0 && i~=1
        drawnow;
        f = f+1;
        saveas(gcf, ['Density_', num2str(bandwidth), '_', num2str(f), '.png']);
        disp([num2str(12*f) ' transcripts are finished.']);
        
        clf;
        set(gcf, 'units', 'normalized', 'position', [.05+.02*f .1-.02*f .8 .8]);
    end
    posDensity = pos(idxName==i,1:2);
    if size(posDensity,1)>2
        [~, density, X, Y]=kde2d_modified...
            (posDensity,2^10,[0 0],...
            floor([imsize(2)/scale/5, imsize(1)/scale/5]),...
            floor([bandwidth/5, bandwidth/5]));
    else
        density = zeros(2^10);
    end
    density = imresize(density, [imsize(1)/5, imsize(2)/5]);
    subplot(3, 4, i-f*12);
    imshow(log2(density), []);
    colormap(gca, parula);
    title([uniName{i}, ' (', num2str(countName(i)), ')']);
end
disp('All transcripts are finished.')
saveas(gcf, ['Density_', num2str(bandwidth), '_', num2str(f+1), '.png']);

