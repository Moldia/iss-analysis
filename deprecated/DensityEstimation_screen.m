%% kernel density estimation of all genes
%  The figure rendering will take quite some time. Be patient!
%  Xiaoyan, 2015-8-10
% deprecated

clear;
close all;

%% parameters
decoded_file = 'input_example\QT_0.45_1e-05_details.csv';
image = 'input_example\HE_10_c1+2+3.jpg';   % important for size
scale = 1;
bandwid = 600; % in original scale

%% transcripts
[name,Pos] = getinsitudata_f(decoded_file);
Pos = correctcoordinates_f(Pos,.2);

% unique transcripts
[name_uni,~,idx_re] = unique(name);
[p,q] = hist(idx_re,unique(idx_re));

%% image size
imgin = imfinfo(image);
Isize = [imgin.Height,imgin.Width];

%% density estimation plot
tic
f = 0;
figure;
set(gcf,'units','normalized','position',[.05+.02*f .1-.02*f .8 .8]);
for i = 1:length(name_uni)
    if mod(i-1,12)==0 && i~=1
        f = f+1;
        disp([num2str(12*f) ' transcripts are finished.']);
        drawnow;
        figure;
        set(gcf,'units','normalized','position',[.05+.02*f .1-.02*f .8 .8]);
    end
    pos_density = Pos(idx_re==i,1:2);
    if size(pos_density,1)>2
        [bandwidth,density,X,Y]=kde2d_X(pos_density,2^10,[0 0],floor([Isize(2)/scale/5 Isize(1)/scale/5]),floor([bandwid/5 bandwid/5]));
    else
        density = zeros(2^10);
    end
    density = imresize(density,[Isize(1)/5 Isize(2)/5]);
    subplot(3,4,i-f*12);
    imshow(density,[]);
    colormap(parula);
    title([name_uni{i} ' (' num2str(p(i)) ')']);
end
disp('All transcripts are finished.')
