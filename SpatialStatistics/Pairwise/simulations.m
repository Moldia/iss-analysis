

imsize = 500;
bandwidth = floor(imsize/15);
fh = fspecial('gaussian',bandwidth*2,bandwidth/5);

[blobs_1,~] = generateABblobs_f(imsize,500,1,20,0);
% blobs_2 = blobs_1;
conv_blobs_1 = imfilter(blobs_1,fh);
conv_blobs_2 = imfilter(blobs_2,fh);
merged_conv = cat(3,conv_blobs_1,conv_blobs_2,conv_blobs_1);

figure;
subplot(121),imshow(merged_conv/max(merged_conv(:)));

% black pixels
step = 0.5;

black = conv_blobs_1==0 & conv_blobs_2==0;
ratio = conv_blobs_1(~black)./conv_blobs_2(~black);
[a,b] = hist(atan(ratio),(0:step:90)/180*pi);

subplot(122),plot(1:length(b),a,'linewidth',2);
xpoint = 1:length(b);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);
