%% Transparent plot for Eph, one color for all
%  transcripts from two or more sections


clear;
close all;

%% parameters
image = 'E:\PROOOJECTS\9_Ephrin\Sample_snapshot\12189_2_8_s1c1.jpg';
scale = 1;

decoded_file_list = {...    % a list of files to read
    'E:\PROOOJECTS\9_Ephrin\Image analysis\12189_A2_EGFR\12189_2_1\Decoded_details.csv',...
    'E:\PROOOJECTS\9_Ephrin\Image analysis\12189_A2_EGFR\12189_8_1\Decoded_details.csv',...
    'E:\PROOOJECTS\9_Ephrin\Image analysis\12189_A2_EGFR\12189_12_1\Decoded_details.csv'};
channel_order_list = {...   % original channel, four channels for each file
    'EGFR mut' 'A2 wt' 'EGFR wt' 'A2 mut';...
    'EGFR mut' 'A2 wt' 'EGFR wt' 'A2 mut';...
    'EGFR mut' 'A2 wt' 'EGFR wt' 'A2 mut'};
name_list = {...
    'EGFR mut';...
    'EGFR mut';...
    'EGFR mut'}; % transcripts to plot, one for each section

% if more than one transcript from one section need to be plotted, pretend
% it to be a transcript from yet another section, i.e. add the same info
% twice to the file list and channel list, add a new name to name list

figure_title = 'EGFR - 12189 all';
show_image = 1;
radius = 100;

%% extract data - different sections
Data = struct;
for s = 1:length(decoded_file_list)
    data = csvread(decoded_file_list{s},1);
    field = ['file' num2str(s)];
    Data.(field).position = data(:,1:2);
    Data.(field).uniname = channel_order_list(s,:);
    Data.(field).uniname_re = data(:,4);
end
    
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
for s = 1:length(decoded_file_list)
    field = ['file' num2str(s)];
    position = Data.(field).position;
    uniname = Data.(field).uniname;
    uniname_re = Data.(field).uniname_re;
    
    idx = find(strcmp(uniname,name_list{s}));
    if isempty(idx)
        warning(['Could not find ' name_list{s} ' in the input file.'])
    else
        x = position(uniname_re==idx,1);
        y = position(uniname_re==idx,2);
        plotx = (repmat(x,1,15)+repmat(circx,length(x),1))';
        ploty = (repmat(y,1,15)+repmat(circy,length(y),1))';
        patch(plotx*scale,ploty*scale,'r','facealpha',.3,'linestyle','none');
    end
    
end
title(figure_title);