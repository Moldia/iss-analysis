%% randomly sample a population of cells
%  cell class/ population id needed
%  based on cell label image
%  spatially map to original cells
%  Xiaoyan, last update 2016-4-5

%%
I = imread('E:\Demo data\NatMeth\SlideA\Stitched\DO_stitched_c1.tif');

load('RenumberedCells.mat');
accensefile = 'accense_output.csv';
data = importdata(accensefile);
data = data.data(:,[1,28]);

%% screening
col = parula(max(data(:,2)));
Ilab = zeros(size(CellLabel,1),size(CellLabel,2),3);
for t = 1:max(data(:,2))
%     if t == 9 || t == 14
        cellclass = t;
        positive = data(data(:,2)==cellclass,1);
        img_positive = double(ismember(CellLabel,positive));
        Ilab = Ilab + cat(3,img_positive*col(t,1),img_positive*col(t,2),img_positive*col(t,3));
%     end
end

figure,
imshow(I*3); hold on;
hi = imshow(Ilab);
set(hi,'alphadata',.5);

% set(gca,'xlim',[2550 3600],'ylim',[180 1000]);

%% subset
Test = [9, 14];

for t = 1:length(Test)
    figure;
    cellclass = Test(t);
    positive = data(data(:,2)==cellclass,1);
    img_positive = uint8(ismember(CellLabel,positive));
    randpos = randperm(length(positive),15);
    subplotpos = [1 2 3 6 7 8 11 12 13 16 17 18 21 22 23];
    
    for i = 1:15
        subplot(5,5,subplotpos(i));
        temp = positive(randpos(i));
        img_positive(CellLabel==temp) = i+1;
        [coordy,coordx] = find(CellLabel==temp);

        imshow(I(max(min(coordy)-30,1):min(max(coordy)+30,size(I,1)),...
            max(min(coordx)-30,1):min(max(coordx)+30,size(I,2)),:),[]);
    end
    subplot(5,5,[4 5 9 10 14 15 19 20 24 25]);
    imshow(I,[]); hold on;
    hi = imshow(label2rgb(img_positive(1:size(I,1),1:size(I,2),:))); 
    set(hi,'alphadata',0.5);
    title(num2str(cellclass));
end


