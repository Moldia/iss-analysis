%% randomly sample marker-positive cells
%  based on cell label image
%  spatially map to original cells
%  plot original signals
%  Xiaoyan, last update 2016-4-5

%%
load('RenumberedCells.mat');
I = imread('wCellBoundaries.tif');
blobfile = 'E:\Demo data\NatMeth\CP_160111_NucleiAlign\Decoding\QT_0.4_1e-20_details.csv';
[names,pos] = getinsitudata_f(blobfile);
[name_uni,~,idx_name] = unique(names);

data = importdata('Cell_blobs.csv');
names = data.textdata(2:end);
count = data.data(2:end,:);

threshold = 2;
Test = {'MUC1' 'KI-67'};

%%
for t = 1:length(Test)
    figure;
    namecell = find(strcmp(names,Test{t}));
    positive = count(count(:,namecell+1)>=threshold,1);
    img_positive = uint8(ismember(CellLabel,positive));
    randpos = randperm(length(positive),6);
    subplotpos = [1 2 3 6 7 8];
    
    nameblob = find(strcmp(name_uni,Test{t}));
    for i = 1:6
        subplot(2,5,subplotpos(i));
        temp = positive(randpos(i));
        img_positive(CellLabel==temp) = i+1;
        [coordy,coordx] = find(CellLabel==temp);

        imshow(I(max(min(coordy)-30,1):min(max(coordy)+30,size(I,1)),...
            max(min(coordx)-30,1):min(max(coordx)+30,size(I,2)),:));
        hold on;
        temp = pos(idx_name==nameblob,:);
        plot(temp(:,1)-max(min(coordx)-30,1),temp(:,2)-max(min(coordy)-30,1)+30,'yo');

    end
    subplot(2,5,[4,5,9,10]);
    imshow(I); hold on;
    hi = imshow(label2rgb(img_positive(1:size(I,1),1:size(I,2),:))); 
    set(hi,'alphadata',0.5);
    title(Test{t});
    
end

