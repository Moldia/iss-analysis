% plot all reads alphabetically
% possible to add transparent image layer
% Xiaoyan, 2017

%% modify here
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
image = 'E:\PROOOJECTS\test_dataset\860502_1_align.png';
scale = .2;     % image scale

%% do not modify

% load and remove NNNN
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);
figure;
plotall(name, pos, image, scale)

% update_legend(gca, unique(name))
add_scalebar()
