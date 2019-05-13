clear;
cd('E:\PROOOJECTS\9_Lung_Multiplex\11692\11692_2_s2_m3');
format compact

%% parameters
Thtop = [.0025 .0018 .001 .0025];
Thback = [.009 .0065 .004 .015];
Rback = [2 2 .8 2];
    
%% prototye
% tissue region
Inuclei = imread('ZENout\11692_2-MIP_s2c1_ORG.tif');
Imask = im2bw(Inuclei,graythresh(Inuclei));
Imask = imdilate(Imask,strel('disk',50));

% load images
I = [];
for c = 1:4
    Itemp = imread(['ZENout\11692_2-MIP_s2c' num2str(c+1) '_ORG.tif']);
    % channel alignment
    if c == 1
        Itemp = [zeros(1,size(Itemp,2));Itemp(1:end-1,2:end),zeros(size(Itemp,1)-1,1)];
    end
    Itemp = double(Itemp)/65525;
    I = cat(3,I,Itemp);
end
disp('image loading finished');

% segmentation
I_thresh = [];
I_threshdil = [];
I_comp = [];
for c = 1:4
    c
    Itemp = I(:,:,c);
    
    % filters
    Iback = imclose(Itemp,strel('disk',10));
    Iback = imdilate(Iback,strel('disk',2));
    Itop = imtophat(Itemp,strel('disk',5));
    
    % thresholding and size filtering
    Ithreshraw = im2bw(Itop,Thtop(c));
    Ithresh = bwareaopen(Ithreshraw,2);
    big = bwareaopen(Ithresh,100);
    Ithresh = Ithresh - big;
    Ithresh(Ithresh<0) = 0;
    
    % signals under autofluorescence
    Icomp = zeros(size(Itemp));
    switch c
        case 1
            Icomp(Itemp>=I(:,:,2)*Rback(c)) = 1;
        case 2
            Icomp(Itemp>=I(:,:,1)*Rback(c)) = 1;
        case 3
            Icomp(Itemp>=I(:,:,1)*Rback(c) & Itemp>=I(:,:,4)*0.5) = 1;
        case 4
            Icomp(Itemp>=I(:,:,1)*Rback(c) & Itemp>=I(:,:,3)*4) = 1;
    end
    Icomp = bwareaopen(Icomp,5);
    
    % strong background structures
    Ithreshback = im2bw(Iback,Thback(c));
    Ithreshback = Ithreshback - Icomp;
    Ithresh(Ithresh<0) = 0;
    Ithreshback = bwareaopen(Ithreshback,150);
    Ithresh = Ithresh - Ithreshback;
    Ithresh(Ithresh<0) = 0;
    
    % remove signals outside tissue region
    Ithresh = Ithresh.*Imask;    
    Ithreshd = imdilate(Ithreshraw,strel('disk',5));   

    % imshow(cat(3,Itemp*50,Ithresh,Ithresh))
    % imshow(cat(3,Itemp*50,Ithreshback,Ithreshback))
    % imshow(cat(3,Itemp*50,Icomp,Icomp))

    I_thresh = cat(3,I_thresh,Ithresh);
    I_threshdil = cat(3,I_threshdil,Ithreshd); 
    I_comp = cat(3,I_comp,Icomp);
end

% signals appearing in multiple channels
Imulti = mean(I_threshdil,3);
Imulti(Imulti>.25) = 1;
Imulti(Imulti==.25) = 0;

% remove signals appearing in multiple channels
count = [];
for c = 1:4
    Itemp = imreconstruct(Imulti,I_threshdil(:,:,c));
    Itemp = imreconstruct(Itemp,I_thresh(:,:,c));
    Itemp = Itemp - I_comp(:,:,c);
    Itemp(Itemp<0) = 0;
    Itemp = I_thresh(:,:,c) - Itemp;
    Itemp(Itemp<0) = 0;
    Itemp = imfill(Itemp,'holes');
    I_thresh(:,:,c) = Itemp;
    
    % count objects
    obj = bwconncomp(Itemp);
    count = [count,obj.NumObjects];
    
end

imshow(cat(3,I_thresh(:,:,2),I_thresh(:,:,1),zeros(size(Itemp))) + ...
    cat(3,I_thresh(:,:,4),zeros(size(Itemp)),I_thresh(:,:,4)) + ...
    cat(3,I_thresh(:,:,3),I_thresh(:,:,3)*.3,zeros(size(Itemp))) + ...
    cat(3,double(Inuclei)/65525,double(Inuclei)/65525,double(Inuclei)/65525))

count

%% read images
I1 = imread('ZENout\11692_2-MIP_s2c2_ORG.tif');
I2 = imread('ZENout\11692_2-MIP_s2c3_ORG.tif');
I3 = imread('ZENout\11692_2-MIP_s2c4_ORG.tif');
I4 = imread('ZENout\11692_2-MIP_s2c5_ORG.tif');

Inuclei = imread('ZENout\11692_2-MIP_s2c1_ORG.tif');


%% scattet plots of multispsectral intensities
% figure;
% subplot(2,2,1),plot(I1(:),I2(:),'.'); axis([0 4000 0 4000])
% subplot(2,2,2),plot(I1(:),I3(:),'.'); axis([0 4000 0 4000])
% subplot(2,2,3),plot(I1(:),I4(:),'.'); axis([0 4000 0 4000])
% subplot(2,2,4),plot(I2(:),I3(:),'.'); axis([0 4000 0 4000])
% 
% figure;
% plot3(I1(:),I2(:),I3(:),'.');
% 
% 
% figure;distributionPlot([double(I2(:))./double(I1(:)),...
%     double(I3(:))./double(I1(:)),...
%     double(I4(:))./double(I1(:)),...
%     double(I3(:))./double(I2(:)),...
%     double(I4(:))./double(I2(:)),...
%     double(I4(:))./double(I3(:))],'histOpt',1.1)
% ylim([0 6]);

%% channel alignment
% I1_red = I1-mean(prctile(I1(:),40));
% % figure,imshow(I1_red,[])
% I2_red = I2-mean(prctile(I2(:),40));
% 
% I  = cat(3,imadjust(I1_red),imadjust(I2_red),zeros(size(I1)));
% figure,imshow(I)

I1 = [zeros(1,size(I1,2));I1(1:end-1,2:end),zeros(size(I1,1)-1,1)];
% imshow(cat(3,I2*100,I1*100,zeros(size(I1))));

%% c2,3,4 to RGB
% I = cat(3,I1*50,I2*50,I3*50);
% 
% imshow(I,[])
% imwrite(I,'merged.tif','tiff','compression','none');
% 
% I_2 = I1+I2+I3;
% imtool(I_2)
% 
% I_3= imtophat(I3,strel('disk',20));
% imtool(I_3)

%% tissue region
Imask = im2bw(Inuclei,graythresh(Inuclei));
Imask = imdilate(Imask,strel('disk',30));
% imshow(Imask);

%% grayscale
I1 = double(I1)/65525;
I2 = double(I2)/65525;
I3 = double(I3)/65525;
I4 = double(I4)/65525;

%% smoothing
% h = fspecial('gaussian',25,3);
% I1smooth = imfilter(I1,h);
% I2smooth = imfilter(I2,h);
% I3smooth = imfilter(I3,h);
% I4smooth = imfilter(I4,h);

%% closing
h = strel('disk',10);
I1close = imclose(I1,h);
I2close = imclose(I2,h);
I3close = imclose(I3,h);
I4close = imclose(I4,h);
disp('closing finished');

h = strel('disk',2);
I1closedil = imdilate(I1close,h);
I2closedil = imdilate(I2close,h);
I3closedil = imdilate(I3close,h);
I4closedil = imdilate(I4close,h);

%% opening
% h = strel('disk',5);
% I1open = imopen(I1,h);
% I2open = imopen(I2,h);
% I3open = imopen(I3,h);
% I4open = imopen(I4,h);

%% tophat
% I1top = I1 - I1open;
% I2top = I2 - I2open;
% I3top = I3 - I3open;
% I4top = I4 - I4open;

h = strel('disk',5);
I1top = imtophat(I1,h);
I2top = imtophat(I2,h);
I3top = imtophat(I3,h);
I4top = imtophat(I4,h);
disp('tophat filtering finished')

%% thresholding - c3
I2thresh = im2bw(I2top,Thtop(2));

% size filtering
% I2threshconn = bwconncomp(I2thresh);
% area = cellfun(@size,I2threshconn.PixelIdxList,'UniformOutput',false);
% area = cell2mat(area');
I2thresh = bwareaopen(I2thresh,2);
big = bwareaopen(I2thresh,100);
I2thresh = I2thresh - big;
I2thresh(I2thresh<0) = 0;

% signals under autofluorescence
I2_comp = zeros(size(I2));
I2_comp(I1(:)*2<=I2(:)) = 1;
I2_comp = bwareaopen(I2_comp,5);   % remove small objects (can be a result of channel misalignment)
% imshow(cat(3,I2*100,I2_comp,I2_comp))

% big background structures
I2threshback = im2bw(I2closedil,Thback(2));
% reconst = imreconstruct(logical(I2thresh),I2threshlow);
% reconst = bwareaopen(reconst,400);
I2threshback = I2threshback - I2_comp;  % remove pixels of hidden signals from the background
I2thresh(I2thresh<0) = 0;
I2threshback = bwareaopen(I2threshback,150);    % background bigger than 150 px
I2thresh = I2thresh - I2threshback;     % remove from the segmented structures
I2thresh(I2thresh<0) = 0;

% outside tissue
I2thresh = I2thresh.*Imask;

% imshow(cat(3,I2*100,I2thresh,I2thresh))
% imshow(cat(3,I2*100,I2threshback,I2threshback))
I2threshd = imdilate(I2thresh,strel('disk',4));

%% thresholding - c2
I1thresh = im2bw(I1top,Thtop(1));

% size filtering
I1thresh = bwareaopen(I1thresh,2);
big = bwareaopen(I1thresh,100);
I1thresh = I1thresh - big;
I1thresh(I1thresh<0) = 0;

% signals under autofluorescence
I1_comp = zeros(size(I1));
I1_comp(I2(:)*1.5<=I1(:)) = 1;
I1_comp = bwareaopen(I1_comp,15);
% imshow(cat(3,I1_comp,I1*50,I1_comp))

% strong background
I1threshback = im2bw(I1closedil,Thback(1));
I1threshback = I1threshback - I1_comp;
I1thresh(I1thresh<0) = 0;
I1threshback = bwareaopen(I1threshback,150);
I1thresh = I1thresh - I1threshback;
I1thresh(I1thresh<0) = 0;

% outside tissue
I1thresh = I1thresh.*Imask;

% imshow(cat(3,I1thresh,I1*100,I1thresh))
% imshow(cat(3,I1threshback,I1*100,I1threshback))
I1threshd = imdilate(I1thresh,strel('disk',4));

%% autofluoresnce & cross-talk
I_comp = I1threshd.*I2threshd;
% imshow(cat(3,I2threshd,I1threshd,Icross))
% imshow(cat(3,I2*100,I1*100,Icross))

%% thresholding - c4
I3thresh = im2bw(I3top,Thtop(3));

% size filtering
I3thresh = bwareaopen(I3thresh,2);
big = bwareaopen(I3thresh,100);
I3thresh = I3thresh - big;
I3thresh(I3thresh<0) = 0;

% signals under autofluorescence
I3_comp = zeros(size(I3));
I3_comp(I1(:)*.8<=I3(:)) = 1;
I3_comp = bwareaopen(I3_comp,5);   % remove small objects (can be a result of channel misalignment)
% imshow(cat(3,I3*100,I3_comp,I3_comp))

% strong background
I3threshback = im2bw(I3closedil,Thback(3));
I3threshback = I3threshback - I3_comp;
I3thresh(I3thresh<0) = 0;
I3threshback = bwareaopen(I3threshback,150);
I3thresh = I3thresh - I3threshback;
I3thresh(I3thresh<0) = 0;

I3thresh = I3thresh.*Imask;

% imshow(cat(3,I3*100,I3thresh,I3thresh))
% imshow(cat(3,I3*100,I3threshback,I3threshback))
I3threshd = imdilate(I3thresh,strel('disk',2));

%% autofluoresnce & cross-talk
I_comp = I1threshd.*I3threshd;
% imshow(cat(3,I1threshd,I3threshd,Icross))
% imshow(cat(3,I3*200,I1*100,Icross))

%% thresholding - c5
I4thresh = im2bw(I4top,Thtop(4));

% size filtering
I4thresh = bwareaopen(I4thresh,2);
big = bwareaopen(I3thresh,100);
I4thresh = I4thresh - big;
I4thresh(I4thresh<0) = 0;

% strong background
I4threshback = im2bw(I4closedil,Thback(4));
I4threshback = bwareaopen(I4threshback,150);
I4thresh = I4thresh - I4threshback;
I4thresh(I4thresh<0) = 0;
I4thresh = I4thresh.*Imask;

% imshow(cat(3,I4*50,I4thresh,I4thresh))
% imshow(cat(3,I4*50,I4threshback,I4threshback))
I4threshd = imdilate(I4thresh,strel('disk',2));

%% autofluoresnce & cross-talk
% Icross = I1threshd.*I4threshd;
% imshow(cat(3,I1threshd,I4threshd,Icross))
% imshow(cat(3,I3*200,I1*100,Icross))

%% remove autofluorescence
Iauto = mean(cat(3,I1threshd,I2threshd,I3threshd,I4threshd),3);
% imshow(Iauto,[.25 1])
Iauto(Iauto>.25) = 1;
Iauto(Iauto==.25) = 0;
% imshow(Iauto)
% imshow(cat(3,I1*100,Iauto,Iauto))

reconst = imreconstruct(Iauto,I1threshd);
reconst = imreconstruct(reconst,I1thresh);
I1thresh = I1thresh - reconst;
I1thresh(I1thresh<0) = 0;

reconst = imreconstruct(Iauto,I2threshd);
reconst = imreconstruct(reconst,I2thresh);
I2thresh = I2thresh - reconst;
I2thresh(I2thresh<0) = 0;

reconst = imreconstruct(Iauto,I3threshd);
reconst = imreconstruct(reconst,I3thresh);
I3thresh = I3thresh - reconst;
I3thresh(I3thresh<0) = 0;

reconst = imreconstruct(Iauto,I4threshd);
reconst = imreconstruct(reconst,I4thresh);
I4thresh = I4thresh - reconst;
I4thresh(I4thresh<0) = 0;

imshow(cat(3,I2thresh,I1thresh,zeros(size(I4thresh))) + ...
    cat(3,I4thresh,zeros(size(I4thresh)),I4thresh) + ...
    cat(3,I3thresh,I3thresh*.3,zeros(size(I4thresh))) + ...
    repmat(double(Inuclei)/65525,1,1,3))


%% tophat - lines kernels to remove autofluorescence
% on original images
% ang = 0:20:180;
% 
% for i = 1:length(ang)-1
%     SE = strel('line',15,ang(i));
%     Itemp = imtophat(I2,SE);
%     if i == 1
%         Itop_line = Itemp;
%     else
%         Itop_line = Itop_line+Itemp;
%     end
% end
% imshow(Itop_line,[]);

%% find autofluorescence - line kernels and reconstruction
% ang = 0:10:180;
% Itest = [];
% Itrack = I2;
% for i = 1:length(ang)
%     SE = strel('line',15,ang(i));
%     Itemp = imerode(I2,SE);   % erosion
%     Itemp = imreconstruct(Itemp,I2open);  % reconstruction
%     Itemp = I2-Itemp; % remove reconstructed item
%     Itest = cat(3,Itest,Itemp);
%     Itrack(Itrack~=Itemp) = 0;
% end
% 
% imshow(Itest(:,:,10:12))
% imtool(mean(Itest,3)*65525)
% 
% Itestmean = Itest - repmat(mean(Itest,3),1,1,4);
% 
% I2rem = I2;
% I2rem(Itrack~=1)=0;
% figure,imshow(I2rem,[])
% 
% Itrack(Itrack~=1)=0;
% 
% Itest = imreconstruct(Itrack,I2rem);

%% LoG
% h = fspecial('log',20,2);
% I1log = imfilter(I1,h);
% I2log = imfilter(I2,h);
% 
% h = fspecial('log',20,4);
% I3log = imfilter(I3,h);
% I4log = imfilter(I4,h);
% 
% %% min and max 
% I1min = min(I1log(:));
% I1max = max(I1log(:));
% 
% I2min = min(I2log(:));
% I2max = max(I2log(:));
% 
% %% invert LoG images
% I1loginv = imcomplement(I1log);
% % I1loginv = I1max-I1log;
% imtool(I1loginv*1E5)
% test = imreconstruct(I1loginv,I1open);
% 
% I2loginv = I2max-I2log;
% imtool(uint16(I2loginv*65525))
% I2thres = im2bw(I2loginv,235/65525);
% I2threslow = im2bw(I2loginv,195/65525);
% I2reconst = imreconstruct(I2thres,I2threslow);

%% autofluorescence segmentation
% Prob_auto = imread('Multispectral_classification\11692_2-MIP_s2c_Probabilities_autofluorescence.png');
% imtool(Prob_auto)
% 
% Mask_auto = im2bw(Prob_auto,247/255);
% figure,imshow(Mask_auto)

%% autoluorescence intensities
% I1_auto = double(I1).*Mask_auto;
% I2_auto = double(I2).*Mask_auto;
% I3_auto = double(I3).*Mask_auto;
% I4_auto = double(I4).*Mask_auto;
% 
% pix_auto = find(Mask_auto(:));
% multispect_auto = cat(2,I1_auto(pix_auto),I2_auto(pix_auto),I3_auto(pix_auto),I4_auto(pix_auto));
% figure,plot(multispect_auto)
% 
% figure;
% subplot(2,2,1),plot(multispect_auto(:,1),multispect_auto(:,2),'.'); axis([0 4000 0 4000])
% subplot(2,2,2),plot(multispect_auto(:,1),multispect_auto(:,3),'.'); axis([0 4000 0 4000])
% subplot(2,2,3),plot(multispect_auto(:,1),multispect_auto(:,4),'.'); axis([0 4000 0 4000])
% subplot(2,2,4),plot(multispect_auto(:,2),multispect_auto(:,3),'.'); axis([0 4000 0 4000])
% 
% figure;
% plot3(multispect_auto(:,1),multispect_auto(:,2),multispect_auto(:,3),'.');
% 
% figure;distributionPlot([multispect_auto(:,2)./multispect_auto(:,1),...
%     multispect_auto(:,3)./multispect_auto(:,1),...
%     multispect_auto(:,4)./multispect_auto(:,1),...
%     multispect_auto(:,3)./multispect_auto(:,2),...
%     multispect_auto(:,4)./multispect_auto(:,2),...
%     multispect_auto(:,4)./multispect_auto(:,3)],'histOpt',1.1)

