
%% parameters
lsize = 200;
bandwidth = 50;
step = 0.25;

% gaussian filter
size = bandwidth*2;
if mod(size,2)==0
    size = size+1;
end
sigma = bandwidth/5;
x = linspace(-size/2,size/2,size);
gaussFilter = exp(-x.^2/(2*sigma^2));
gaussFilter = gaussFilter/sum(gaussFilter);

x = 1:lsize;

%% with specific distance, blobs 1 fixed
mean_blob_1 = 80;

figure;
blobs_1 = zeros(lsize,1);
blobs_1(mean_blob_1) = 1;
conv_blobs_1 = conv(blobs_1,gaussFilter,'same');
gaussian_analog = exp((-(x-mean_blob_1).^2)/(2*sigma^2));
gaussian_analog = gaussian_analog/sum(gaussian_analog);

subplot(3,4,1:2);
plot(conv_blobs_1/max(gaussFilter),'linewidth',2);
hold on;
plot(gaussian_analog/max(gaussFilter),'--','linewidth',2);
box off

for dist = 2:20:100
    
    blobs_2 = zeros(lsize,1);
    blobs_2(mean_blob_1+dist) = 1;
    conv_blobs_2 = conv(blobs_2,gaussFilter,'same');
    
    subplot(3,4,3:4);
    hold on;
    plot(conv_blobs_2/max(gaussFilter),'linewidth',2);
    
    % color ratio
    black = conv_blobs_1==0 & conv_blobs_2==0;
    ratio = conv_blobs_2(~black)./conv_blobs_1(~black);
    [a,b] = hist(atan(ratio),(0:step:90)/180*pi);
    
    subplot(3,1,2)
    hold on;
    plot(1:length(b),a,'linewidth',2);
   
%     gaussian_analog = gaussian_analog/sum(gaussian_analog);
    ratio_analog = exp((-((0:0.0001:lsize)-(mean_blob_1+dist)).^2+((0:0.0001:lsize)-mean_blob_1).^2)/(2*sigma^2));
    [c,~] = hist(atan(ratio_analog),(0:step:90)/180*pi);
    subplot(3,1,3)
    hold on;
    plot(1:length(b),c,'linewidth',2);

end

xpoint = 1:length(b);

subplot(3,4,3:4);
plot(gaussian_analog/max(gaussFilter),'--','linewidth',2);
subplot(3,1,2);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);
legend(dist);
subplot(3,1,3);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);
legend(dist);

%% with specific distance, pair
figure;
col = parula(4);

i = 0;
for dist = 22:2:26
    i = i+1;
    mean_blob_1 = 100-dist/2;
    blobs_1 = zeros(lsize,1);
    blobs_1(mean_blob_1) = 1;
    conv_blobs_1 = conv(blobs_1,gaussFilter,'same');
    gaussian_analog = exp((-(x-mean_blob_1).^2)/(2*sigma^2));
    gaussian_analog = gaussian_analog/sum(gaussian_analog);
    
    subplot(3,1,1);
    hold on;
    plot(gaussian_analog/max(gaussFilter),'--','linewidth',2,'color',col(i,:));
    box off
    
    blobs_2 = zeros(lsize,1);
    blobs_2(mean_blob_1+dist) = 1;
    conv_blobs_2 = conv(blobs_2,gaussFilter,'same');
    
    subplot(3,1,1);
    hold on;
    plot(conv_blobs_2/max(gaussFilter),'linewidth',2,'color',col(i,:));
    
    % color ratio
    black = conv_blobs_1==0 & conv_blobs_2==0;
    ratio = conv_blobs_2(~black)./conv_blobs_1(~black);
    [a,b] = hist(atan(ratio),(0:step:90)/180*pi);
    
    subplot(3,1,2)
    hold on;
    plot(1:length(b),a,'linewidth',2,'color',col(i,:));
   
%     gaussian_analog = gaussian_analog/sum(gaussian_analog);
    ratio_analog = exp((-((0:0.0001:lsize)-(mean_blob_1+dist)).^2+((0:0.0001:lsize)-mean_blob_1).^2)/(2*sigma^2));
    [c,~] = hist(atan(ratio_analog),(0:step:90)/180*pi);
    subplot(3,1,3)
    hold on;
    plot(1:length(b),c,'linewidth',2,'color',col(i,:));

end

xpoint = 1:length(b);

subplot(3,1,2);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);
subplot(3,1,3);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);

%% with specific distance, with cross-talk, blobs 1 fixed
mean_blob_1 = 80;

figure;
blobs_1 = zeros(lsize,1);
blobs_1(mean_blob_1) = 1;
conv_blobs_1 = conv(blobs_1,gaussFilter,'same');
gaussian_analog = exp((-(x-mean_blob_1).^2)/(2*sigma^2));
gaussian_analog = gaussian_analog/sum(gaussian_analog);

subplot(3,4,1:2);
plot(conv_blobs_1/max(gaussFilter),'linewidth',2);
hold on;
plot(gaussian_analog/max(gaussFilter),'--','linewidth',2);
ylim([0 2]);
box off

for dist = 2:20:100
    
    blobs_2 = zeros(lsize,1);
    blobs_2(mean_blob_1+dist) = 1;
    blobs_2 = blobs_2 + blobs_1;
    conv_blobs_2 = conv(blobs_2,gaussFilter,'same');
    
    subplot(3,4,3:4);
    hold on;
    plot(conv_blobs_2/max(gaussFilter),'linewidth',2);
    
    % color ratio
    black = conv_blobs_1==0 & conv_blobs_2==0;
    ratio = conv_blobs_2(~black)./conv_blobs_1(~black);
    [a,b] = hist(atan(ratio),(0:step:90)/180*pi);

    subplot(3,1,2)
    hold on;
    plot(1:length(b),a/nnz(~black),'linewidth',2);
    
    ratio_analog = exp((-((0:0.0001:lsize)-(mean_blob_1+dist)).^2+((0:0.0001:lsize)-mean_blob_1).^2)/(2*sigma^2));
    ratio_analog = ratio_analog+1;
    [c,~] = hist(atan(ratio_analog),(0:step:90)/180*pi);
    subplot(3,1,3)
    hold on;
    plot(1:length(b),c,'linewidth',2);
    
end

xpoint = 1:length(b);

subplot(3,4,3:4);
plot(gaussian_analog/max(gaussFilter),'--','linewidth',2);
subplot(3,1,2);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);
legend(dist);
subplot(3,1,3);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);
legend(dist);

%% random blobs
num_blos = 20;

% generate random distribution of siganls
blobs_1 = randi(lsize,num_blos,1);
blobs_1 = accumarray(blobs_1,1,[lsize,1]);
conv_blobs_1 = conv(blobs_1,gaussFilter,'same');

blobs_2 = randi(lsize,num_blos,1);
blobs_2 = accumarray(blobs_2,1,[lsize,1]);
conv_blobs_2 = conv(blobs_2,gaussFilter,'same');

subplot(2,1,1),plot(conv_blobs_1/max(gaussFilter));
hold on;
plot(conv_blobs_2/max(gaussFilter));
stem(blobs_1,'.');
stem(blobs_2,'.');

black = conv_blobs_1==0 & conv_blobs_2==0;
ratio = conv_blobs_2(~black)./conv_blobs_1(~black);
[a,b] = hist(atan(ratio),(0:step:90)/180*pi);
xpoint = 1:length(b);
subplot(2,1,2)
hold on;
plot(1:length(b),a/nnz(~black),'linewidth',2);
set(gca,'XLim',[0 length(b)+1],'XTick',linspace(1,xpoint(end),7),'XTickLabel',0:15:90);

figure;
plot(conv_blobs_1,conv_blobs_2,'.')
