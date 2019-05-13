%% Transparent plot for Eph, one color for all
%  Transcripts from the same section

%% parameters
decoded_file = 'E:\PROOOJECTS\9_Ephrin\Image analysis\8128_A2_EGFR\8128_3_1\Decoded_details.csv';
image = 'E:\PROOOJECTS\9_Ephrin\Sample_snapshot\8128_3_5_s1c1.jpg';
scale = 1;
channel_order = {'EGFR mut' 'A2 wt'	'EGFR wt' 'A2 mut'};    % original channel order
names = {'EGFR mut' 'EGFR wt'}; % transcripts to plot
figure_title = 'EGFR';
show_image = 1;
radius = 100;

%% extract data
data = csvread(decoded_file,1);
position = data(:,1:2);

% unique transcripts
uniname = channel_order;
uniname_re = data(:,4);

%% circles
alph = linspace(0,2*pi,15);
circx = radius*cos(alph);
circy = radius*sin(alph);

%% prepare image
figure;
if show_image
    img = imread(image);
    imshow(img);
else
    set(gca,'YDir','reverse','XTick',[],'YTick',[]);
    axis image;
    axis off;
end

%% plot
hold on;
for i = 1:length(names)
    idx = find(strcmp(uniname,names{i}));
    if isempty(idx)
        warning(['Could not find ' names{i} ' in the input file.'])
    else
        x = position(uniname_re==idx,1);
        y = position(uniname_re==idx,2);
        plotx = (repmat(x,1,15)+repmat(circx,length(x),1))';
        ploty = (repmat(y,1,15)+repmat(circy,length(y),1))';
        patch(plotx*scale,ploty*scale,'r','facealpha',.3,'linestyle','none');
    end
end
title(figure_title);