function Plotting_tile_Sequencing(output_directory_decode,...
    output_filename_afterQT_prefix,...
    single_tile_YN,tile_background_image,...
    multiple_tiles_YN,background_image_tile_prefix,background_image_tile_suffix,...
    continuous_tile_id_YN,max_tile_number,subset_of_tiles_YN,tile_id_to_plot,...
    taglist_plot_tile,tile_use_old_style_taglist_YN,symbol_size_tile,...
    output_directory_plot_tile,output_filename_plot_tile_prefix,...
    plot_tile_exclude_NNNN_YN,plot_tile_reads_beforeQT_YN,...
    plot_tile_based_on_group_YN,tile_use_default_symbol_list,...
    plot_tile_base1_general_stain,varargin)

% Plot reads based on the given taglist_plot_tile.
%
% The background images specified has to be full size, but can be of lower
% quality.
%
% Plotting_tile_Sequencing v2.2.1
% Xiaoyan, 2015-6-20

%% initiate
disp('Initiating Plotting_tile_Sequencing.');
output_directory_decode = [output_directory_decode '\'];
load([output_directory_decode output_filename_afterQT_prefix '.mat']);

output_directory_plot_tile = [output_directory_plot_tile '\'];


%% check tiled images
disp('checking images..');
if single_tile_YN
    if exist(tile_background_image,'file')==2
    else
        error(['Could not find file ' tile_background_image '.']);
    end
elseif multiple_tiles_YN
    if continuous_tile_id_YN
        for i = 1:max_tile_number
            image_name = [background_image_tile_prefix ...
                        num2str(i) background_image_tile_suffix];
            if exist(image_name,'file')==2
            else
                error(['Could not find file ' image_name '.']);
            end
        end
    elseif subset_of_tiles_YN
        for i = 1:length(tile_id_to_plot)
            image_name = [background_image_tile_prefix...
                num2str(tile_id_to_plot(i))...
                background_image_tile_suffix];
            if exist(image_name,'file')==2
            else
                error(['Could not find file ' image_name '.']);
            end
        end
    else
        error('Choose at least one type of multiple tile plotting.');
    end
else
    error('Choose at least one type of tile plotting.');
end
                

%% check and make taglist
[taglist,expected_list,~,~,exp_tags,exp_groups] =  ...
    maketag(taglist_plot_tile,num_hybs,cycle5_empty_YN,tile_use_old_style_taglist_YN,...
    plot_tile_based_on_group_YN,abnormal_sequencing_YN,sequencing_order);

%% prepare symbol list
if tile_use_default_symbol_list
    sym = default_symbolist;
else
    n = size(taglist,2);
    if n==3
        sym = taglist(:,3);
    else
        disp('No symbol column detected in taglist_plot. Continue with default symbols.');
        sym = default_symbolist;
    end
end

%% prepare the read list
if plot_tile_reads_beforeQT_YN
    list = allbt; tile_x = tile_x_allbt; tile_y = tile_y_allbt;
    list_tile = tile_allbt;
else
    list = allqt; tile_x = tile_x_allqt; tile_y = tile_y_allqt;
    list_tile = tile_allqt;
end

%% merge reads with same gene name/gene group
[list,expected_list,unitag] = ...
    groupreadsforplotting(list,expected_list,exp_tags,exp_groups,plot_tile_based_on_group_YN,sym);

%% open images and start plotting
disp('start plotting..');

n=0; data =[]; plot_tag=[];
if single_tile_YN
    if length(unique(list_tile))==1
        I = imread(tile_background_image);
        f=figure; imshow(I,[0 2000]); hold on; axis equal;
        colormap(gray);
        
        if plot_tile_base1_general_stain
            plot(tile_x,tile_y,sym{i},'markersize',symbol_size_tile);
            count = length(list);
            n = n+length(list);
            plot_tag = 'general stain';
        else
            A=[]; B=[]; count=[];
            for i = 1:length(expected_list)
                a = find(list == expected_list(i));
                if a
                    plot(tile_x(a),tile_y(a),sym{i},'markersize',symbol_size_tile);
                    n = n+length(a);
                    plot_tag = [plot_tag unitag(i)];
                    A = [A; a]; count = [count; length(a)];
                else
                    plot_tag_empty = [plot_tag_empty unitag(i)];
                    B = [B; 0];
                end
            end
            if plot_tile_exclude_NNNN_YN
            else
                if length(A) ~= length(list)
                    tile_x_N = removerows(tile_x,'ind',A);
                    tile_y_N = removerows(tile_y,'ind',A);
                    plot(tile_x_N,tile_y_N,'ch');
                    plot_tag = [plot_tag 'NNNN'];
                    count = [count; length(tile_x_N)];
                    n = n+length(tile_x_N);
                else
                    plot_tag_empty = [plot_tag_empty 'NNNN'];
                    B = [B; 0];
                end
            end
        end  
        
        if isempty(B)
            data = count; rowname = [plot_tag];
        else
            data = [count;B]; rowname = [plot_tag, plot_tag_empty];
        end
        
        set(gca,'YDir','reverse');axis off;
        set(gca,'unit','normalized','position',[0.05 0.05 0.9 0.9]);
        set(gcf,'Visible','on','units','normalized');
        set(f,'name','single tile plotting');
        legend(plot_tag);
        
        if exist(output_directory_plot_tile, 'dir')  
        else mkdir (output_directory_plot_tile);
        end 
        imgout = [output_directory_plot_tile ...
            output_filename_plot_tile_prefix...
            '_SingleTile'];
        saveas(gcf, [imgout '.png'], 'png');
        saveas(gcf, [imgout '.fig'], 'fig');
   
    else
        error(['More than one tile IDs are found in the decoding file. '...
            'Use other plotting mode than single tile plotting']);
    end
    
    
elseif multiple_tiles_YN
    
    if continuous_tile_id_YN 
        for i = 1:max_tile_number
            list_sub = list(list_tile==i);
            x_sub = tile_x(list_tile==i);
            y_sub = tile_y(list_tile==i);
            if length(list_sub)
                image_name = [background_image_tile_prefix...
                    num2str(i) background_image_tile_suffix];
                I = imread(image_name);
                f=figure; imshow(I,[0 2000]); hold on; axis equal;
                colormap(gray);
                
                if plot_tile_base1_general_stain
                    plot(x_sub,y_sub,sym{i},'markersize',symbol_size_tile);
                    count = length(x_sub);
                    n = n+length(x_sub);
                    plot_tag = 'general stain';
                else
                    A=[]; B=[]; count=[];
                    for j = 1:length(expected_list)
                        a = find(list_sub == expected_list(j));
                        if a
                            plot(x_sub(a),y_sub(a),sym{j},'markersize',symbol_size_tile);
                            n = n+length(a);
                            plot_tag = [plot_tag unitag(j)];
                            A = [A; a];
                        end
                        count = [count; length(a)];
                    end
                    
                    if plot_tile_exclude_NNNN_YN
                    else
                        if length(A) ~= length(list)
                            x_sub_N = removerows(x_sub,'ind',A);
                            y_sub_N = removerows(y_sub,'ind',A);
                            plot(x_sub_N,y_sub_N,'ch');
                            plot_tag = [plot_tag 'NNNN'];
                            n = n+length(x_sub_N);
                        end
                        count = [count; length(x_sub_N)];
                    end
                    
                end
                
                set(gca,'YDir','reverse');axis off;
                set(gca,'unit','normalized','position',[0.05 0.05 0.9 0.9]);
                set(gcf,'Visible','on','units','normalized');
                legend(plot_tag);
                set(f,'name',['tile ' num2str(i) ' plotting']);
                data = [data,count];
                
                if exist(output_directory_plot_tile, 'dir')
                else mkdir (output_directory_plot_tile);
                end 
                imgout = [output_directory_plot_tile ...
                    output_filename_plot_tile_prefix...
                    '_Tile' num2str(i)];
                saveas(gcf, [imgout '.png'], 'png');
                saveas(gcf, [imgout '.fig'], 'fig');
                
            else
                if plot_tile_base1_general_stain
                    data = [data,0];
                else
                    if plot_tile_exclude_NNNN_YN
                        data = [data,zeros(length(expected_list)+1,1)];
                    else
                        data = [data,zeros(length(expected_list),1)];
                    end
                end
            end
        end
        
    elseif subset_of_tiles_YN
        for i = 1:length(tile_id_to_plot)
            list_sub = list(list_tile==tile_id_to_plot(i));
            x_sub = tile_x(list_tile==tile_id_to_plot(i));
            y_sub = tile_y(list_tile==tile_id_to_plot(i));
            if length(list_sub)
                image_name = [background_image_tile_prefix...
                    num2str(tile_id_to_plot(i))...
                    background_image_tile_suffix];
                I = imread(image_name);
                f=figure; hold on; imagesc(I); axis equal;
                colormap(gray);
                
                if plot_tile_base1_general_stain
                    plot(x_sub,y_sub,sym{i},'markersize',symbol_size_tile);
                    count = length(x_sub);
                    n = n+length(x_sub);
                    plot_tag = 'general stain';
                else
                    A=[]; B=[]; count=[];
                    for j = 1:length(expected_list)
                        a = find(list_sub == expected_list(j));
                        if a
                            plot(x_sub(a),y_sub(a),sym{j},'markersize',symbol_size_tile);
                            n = n+length(a);
                            plot_tag = [plot_tag unitag(j)];
                            A = [A; a];
                        end
                        count = [count; length(a)];
                    end
                    if plot_tile_exclude_NNNN_YN
                    else
                        if length(A) ~= length(list)
                            x_sub_N = removerows(x_sub,'ind',A);
                            y_sub_N = removerows(y_sub,'ind',A);
                            plot(x_sub_N,y_sub_N,'ch');
                            plot_tag = [plot_tag 'NNNN'];
                            n = n+length(tile_x_N);
                        end
                        count = [count; length(tile_x_N)];
                    end
                end
                
                set(gca,'YDir','reverse');axis off;
                set(gca,'unit','normalized','position',[0.05 0.05 0.9 0.9]);
                set(gcf,'Visible','on','units','normalized');
                legend(plot_tag);
                set(f,'name',['tile ' num2str(tile_id_to_plot(i)) ' plotting']);
                data = [data,count];
                
                if exist(output_directory_plot_tile, 'dir')  
                else mkdir (output_directory_plot_tile);
                end 
                imgout = [output_directory_plot_tile ...
                    output_filename_plot_tile_prefix...
                    '_Tile' num2str(tile_id_to_plot(i))];
                saveas(gcf, [imgout '.png'], 'png');
                saveas(gcf, [imgout '.fig'], 'fig');

            else
                if plot_tile_base1_general_stain
                    data = [data,0];
                else
                    if plot_tile_exclude_NNNN_YN
                        data = [data,zeros(length(expected_list)+1,1)];
                    else
                        data = [data,zeros(length(expected_list),1)];
                    end
                end
            end
        end
        
    end
    
end



% table to show counts
%====================================================================
if single_tile_YN
    rown = plot_tag;
    columnn = 'count';
elseif multiple_tiles_YN
    rown = []; columnn = [];
    if continuous_tile_id_YN
        for i = 1:max_tile_number
            columnn = [columnn; ['tile' num2str(i)]];
        end
        if plot_tile_base1_general_stain
            rown = 'general_stain';
        else
            if plot_tile_exclude_NNNN_YN
                rown = unitag;
            else
                rown = [unitag, 'NNNN'];
            end
        end
    else
        for i = 1:length(tile_id_to_plot)
            columnn = [columnn; ['tile' num2str(tile_id_to_plot(i))]];
        end
        if plot_tile_base1_general_stain
            rown = 'general_stain';
        else
            if plot_tile_exclude_NNNN_YN
                rown = unitag;
            else
                rown = [unitag, 'NNNN'];
            end
        end
    end
end
        
f2=figure;
set(f2,'units','normalized','position',[0.3 0.2 0.2 0.6]);
h = uitable(f2,'data',data,'ColumnName',columnn,'RowName',rown);
set(h,'units','normalized','position',[0 0 1 1]);

show = ['Number of blobs plotted: ' num2str(n)];
disp(show);
fprintf('Plotting finished.\n\n'); 

end


