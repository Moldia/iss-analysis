function plotall_f(imgfile, scale, blobfile)

try
    I = imread(imgfile);    
    figure;
    imshow(I,[]);
    hold on;
catch ME
    if strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
        figure;
        hold on;
        axis image;
        set(gca, 'YDir', 'reverse');
        axis(gca, 'off');
    end
end

[names, pos] = getinsitudata_f(blobfile);
[name_uni, ~, idx_name] = unique(names);
pos = correctcoordinates_f(pos, scale);

% remove NNNN
idx_NNNN = find(strcmp(name_uni,'NNNN'));
if ~isempty(idx_NNNN)
    names = names(idx_name ~= idx_NNNN);
    pos = pos(idx_name ~= idx_NNNN, :);
    [name_uni, ~, idx_name] = unique(names);
end

sym = repmat(default_symbolist, ceil(length(name_uni)/length(default_symbolist)), 1);

for i = 1:length(name_uni)
    plot(pos(idx_name==i,1), pos(idx_name==i,2), sym{i});
end

legend(name_uni,'location','NorthEastOutside','color',[.6 .6 .6]);


end
