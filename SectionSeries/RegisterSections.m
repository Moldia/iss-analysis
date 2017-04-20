%% align a series of floating image on top of a reference image
%  need to give rotation angle and translation matrix

clear;
close all;
%%
mkdir('Registered_images_20%');
mkdir('Registered_positions');

%% Import the data
[~, ~, raw] = xlsread('E:\PROOOJECTS\9_Ephrin\Sample_overview.xlsx','registration');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,2:4);
raw = raw(:,[1,5:10]);

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
sample = data(2:end,1);
section = cellVectors(2:end,1);
image = cellVectors(2:end,2);
mutationgroup = cellVectors(2:end,3);
scaling = data(2:end,2);
rotation = data(2:end,3);
xbefore = data(2:end,4);
ybefore = data(2:end,5);
xafter = data(2:end,6);
yafter = data(2:end,7);

clear R raw data cellVectors


%% organize data
[sample_uni,~,sample_re] = unique(sample);
scene = cellfun(@(image) strsplit(image,'/'),image,'UniformOutput',false);
xleft = xbefore - xafter;
yup = ybefore - yafter;

%% register and save images
for i = 1:length(sample_uni)
    % registration info about sample x
    sub_reg = [scaling(sample_re==i),rotation(sample_re==i),xleft(sample_re==i),yup(sample_re==i)];
    sub_section = [section(sample_re==i,:),image(sample_re==i,:),mutationgroup(sample_re==i)];
    sub_scene = scene(sample_re==i);
    
    % remove HE and TP53_KRAS sections
    sub_reg = sub_reg(~(strcmp(sub_section(:,3),'HE') | strcmp(sub_section(:,3),'TP53_KRAS')),:);
    sub_scene = sub_scene(~(strcmp(sub_section(:,3),'HE') | strcmp(sub_section(:,3),'TP53_KRAS')),:);
    sub_section = sub_section(~(strcmp(sub_section(:,3),'HE') | strcmp(sub_section(:,3),'TP53_KRAS')),:);

    % find reference image
    refidx = sub_reg(:,1)==1 & sub_reg(:,2)==0 & sub_reg(:,3)==0 & sub_reg(:,4)==0;
    
    if nnz(refidx)
        refdir = ['F:\EPH_project\Image analysis\'...
            num2str(sample_uni(i)) '_' sub_section{refidx,3} '\'...
            sub_section{refidx,1} '\ZENout'];
        filels = ls(refdir);
        refimg = [refdir '\' filels(3,:)];
        ref = imread(refimg);
        ref = imresize(ref,0.2);
        size_ref = size(ref);
        imwrite(ref,['Registered_images_20%\' sub_section{refidx,1} '_ref.png'],'png');
        
        copyfile(['E:\PROOOJECTS\9_Ephrin\Image analysis\'...
            num2str(sample_uni(i)) '_' sub_section{refidx,3} '\'...
            sub_section{refidx,1} '\Decoded_details.csv'],...
            ['Registered_positions\' sub_section{refidx,1} '_registered.csv']);
        
        mergecheck = ref;
        
        for j = 1:size(sub_reg,1)
            if j ~= find(refidx) && ~isnan(sub_reg(j,2))
                flodir = ['F:\EPH_project\Image analysis\'...
                    num2str(sample_uni(i)) '_' sub_section{j,3} '\'...
                    sub_section{j,1} '\ZENout'];
                filels = ls(flodir);
                floimg = [flodir '\' filels(3,:)];
                flo = imread(floimg);
                flo = imresize(flo,0.2*sub_reg(j,1));
                size_original = size(flo);
                flo = imrotate(flo,-sub_reg(j,2));
                size_rotate = size(flo);
                
                xextra = size_rotate(2)-size_original(2)+size_original(2)-size_ref(2);
                yextra = size_rotate(1)-size_original(1)+size_original(1)-size_ref(1);
                
                % position
                posfile = ['E:\PROOOJECTS\9_Ephrin\Image analysis\'...
                    num2str(sample_uni(i)) '_' sub_section{j,3} '\'...
                    sub_section{j,1} '\ImageAnalysis.mat'];
                % transform positions
                new_pos = transformpositions_f(posfile,sub_reg(j,1),-sub_reg(j,2),....
                    round(sub_reg(j,4)/5+yextra/2,0),...
                    round(sub_reg(j,3)/5+xextra/2,0),...
                    size(flo),5);
                
                flo = transformimage_f_2(flo,['Registered_images_20%\' sub_section{refidx,1} '_ref.png'],...
                    round(sub_reg(j,4)/5+yextra/2,0),...
                    round(sub_reg(j,3)/5+xextra/2,0),...
                    ['Registered_images_20%\' sub_section{j,1} '_registered.png'],new_pos,0.2);
                
                fid = fopen(['Registered_positions\' sub_section{j,1} '_registered.csv'],'w');
                fprintf(fid,'x_pos,y_pos,size,channel\n');
                fprintf(fid,'%d,%d,%d,%d\n',new_pos');
                fclose(fid);
                
                mergecheck = mergecheck + flo;
                
            end
        end
        h = figure;
        imshow(mergecheck,[]);
        waitfor(h)

    else
        error([sample_uni{i} ' cannot find the reference image']);
    end
end
    