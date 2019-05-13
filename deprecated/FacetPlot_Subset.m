%% plot each transcript separately
%  Xiaoyan, 2015-1-10

clear;
close all;
drawnow;

%% parameters
decoded_file = 'input_example\QT_0.45_1e-05_details.csv';
output_folder = ''; 
transcript_to_plot = {'mACTB' 'mCcr4' 'mCcr6'};

%% transcripts
[name,Pos] = getinsitudata_f(decoded_file);

% unique transcripts
[name_uni,~,idx_re] = unique(name);
[p,q] = hist(idx_re,unique(idx_re));

%% plotting
f = 0;
for i = 1:length(transcript_to_plot)
    if mod(i-1,12)==0
        figure;
        f = f+1;
        set(gcf,'units','normalized','position',[.05+.02*f .1-.02*f .8 .8],...
            'visible','off');
        s = 0;
    end
    
    s = s+1;
    subplot(3,4,s);
    idx_plot = find(strcmp(name_uni,transcript_to_plot{i}));
    if isempty(idx_plot)
        title([transcript_to_plot{i} ' (' num2str(0) ')']);
    else
        plot(Pos(idx_re==idx_plot,1),Pos(idx_re==idx_plot,2),'.');
        set(gca,'YDir','reverse');
        axis equal; axis off;
        title([transcript_to_plot{i} ' (' num2str(p(idx_plot)) ')']);
    end
end

%% output
if ~exist(output_folder,'dir')
    mkdir(output_folder);
end

if ~isempty(output_folder)
    if  output_folder(end)~='\'
        output_folder = [output_folder '\'];
    end
end

disp('saving images..');
while f>=1
    figure(f);
    set(gcf,'visible','on');
    saveas(gcf,[output_folder 'FacetPlot_Subset_' num2str(f) '.png'],'png');
    f = f-1;
end