% ISS image analysis workshop, 170614
% Xiaoyan
% tested on MATLAB R2016b

close all; clear;


%% FISH-quant
load ICP/Spots
I = double(imread('ICP/anchor_image.tif'));
isolated = find(GoodIsolated);
averageSpot = zeros(13);

for i = 1:length(isolated)
    spotpos = round(GoodGlobalYX(isolated(i),:));
    try
        spot = I(spotpos(1)-6:spotpos(1)+6, spotpos(2)-6:spotpos(2)+6);
    end
    if nnz(averageSpot)==0
        averageSpot = spot;
    else
        averageSpot = mean(cat(3, averageSpot, spot), 3);
    end

end

figure; bar3(averageSpot);

%% pixel-calling
clear; close all;
load WS170614 Itop

figure(1); clf;
set(gcf, 'name', 'Aligned sequencing images', 'units', 'normalized', 'position', [0 0 1 1]); 
channels = {'DAPI', 'anchor', 'T', 'G', 'C' 'A'};

Ax = [];
for b = 1:4
     for c = 2:6
         ax = subplot(4,6,(b-1)*6+c); Ax = [Ax, ax];
         imshow(Itop{b,c}, []);
         title(['base' num2str(b) ' ' channels{c}]);
     end
end
linkaxes(Ax, 'xy');
pause()

% pixel calling
Ipixelcall = cell(4,1);
for b = 1:4
    [~, I] = max(cat(3, Itop{b,3:6}), [], 3);
    ax = subplot(4,6,(b-1)*6+1); Ax = [Ax, ax];
    imagesc(I);     
    axis image; axis off; set(gca, 'ydir', 'reverse');
    colormap lines(5);
    Ipixelcall{b} = I;
end
linkaxes(Ax, 'xy');
pause()

% remove noise
for b = 1:4
    I1 = Ipixelcall{b};
    
    Ibw = im2bw(Itop{1,2}, .0025);
    Ibw = padimg(Ibw, size(I1,2)-size(Ibw,2), size(I1,1)-size(Ibw,2));

    I2 = zeros(size(I1));
    for c = 1:4
        I = I1.*double(I1==c).*double(Ibw);
%         I = I1.*double(I1==c);
        cc = bwconncomp(logical(I));
        ccsz = cellfun(@length, cc.PixelIdxList);
        cc = cc.PixelIdxList(ccsz>4 & ccsz<200);
        
        I = false(size(I));
        I(cat(1, cc{:})) = true;
        I = bwlabel(I, 8);
        props = regionprops(I, 'MinorAxisLength');
        props = cat(1, props.MinorAxisLength);
        I(ismember(I, find(props<2))) = 0;

        I2(I~=0) = c;
    end
    
    subplot(4,6,(b-1)*6+1);
    imagesc(I2); 
    axis image; axis off; set(gca, 'ydir', 'reverse');
    colormap lines(5);
end
        
