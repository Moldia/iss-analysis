%% connect objects cut by tiling lines and assign blob to its parent cell
%  based on segmented object images (label)
%  Xiaoyan, last update 161106
% deprecated


% parpool
tilepos = getposition('..\ZENout\Tiled_161021_Kv21.csv');

tile_size = tilepos(2,2)-tilepos(1,2);
tiles_x = range(tilepos(:,2))/tile_size + 1;
tiles_y = range(tilepos(:,3))/tile_size + 1;

digit_num = 3;

mkdir('Stitched\');
mkdir('ParentCell\')

%% read files and relabel nuclei objects (parallel)
Ilabeled = cell(tiles_y,tiles_x);
toremove = [];

for i = 1:tiles_y
    parfor j = 1:tiles_x
        tilenum = tiles_x*(i-1)+j;
        if ismember(tilenum,toremove)
            I = zeros(tile_size,tile_size,'uint16');
        else
            tile = num2str(tilenum);
            tilenum = floor(log10(tilenum))+1;
            while tilenum < digit_num
                tile = ['0', tile];
                tilenum = tilenum+1;
            end
            I = imread(['ObjectImages\Kv21_' tile '.tif']);
            uni_I = unique(I(:));
            label = 0;
            I = I(:);
            if length(uni_I) == 1
                I(:) = 0;
            else
                for k = 1:length(uni_I)
                    I(I==uni_I(k)) = label;
                    label = label+1;
                end
            end
            
            I = reshape(I,tile_size,[]);
        end
        Ilabeled(i,j) = {I};
        
    end
    
end

%% stitch relabeled images together
Istitched = zeros(tiles_y*tile_size,tiles_x*tile_size,'uint32');
Label = 0;
for i = 1:tiles_y
    for j = 1:tiles_x
        I = Ilabeled{i,j};
        I = uint32(I);
        I = I(:);
        
        I(I~=0) = I(I~=0) + Label;
        
        I = reshape(I,tile_size,[]);
        Istitched((i-1)*tile_size+1:i*tile_size,(j-1)*tile_size+1:j*tile_size) = I;
        
        if max(I(:))
            Label = max(I(:));
        end
        
    end
end
clear Ilabeled I Label

%% find objects cut by tiling lines and correct
match_x = [];
for i = 1:tiles_x-1
    temp = Istitched(:,tile_size*i:tile_size*i+1);
    temp = temp(temp(:,1)~=0 & temp(:,2)~=0,:);
    temp = unique(temp,'rows');
    match_x = [match_x;temp];
    temp1 = unique(temp(:,1));
    for j = 1:length(temp1)
        if nnz(temp(:,1)==temp1(j)) > 1
            temp2 = temp(temp(:,1)==temp1(j),2);
            temp2 = unique(temp2);
            match_x =[match_x; [temp2(1:end-1),repmat(temp2(end),length(temp2)-1,1)]];
        end
    end
end

for i = 1:size(match_x,1)
    i
    Istitched(Istitched==match_x(i,2)) = match_x(i,1);
end

match_y = [];
for i = 1:tiles_y-1
    temp = Istitched(tile_size*i:tile_size*i+1,:);
    temp = temp(:,temp(1,:)~=0 & temp(2,:)~=0);
    temp = temp';
    temp = unique(temp,'rows');
    match_y = [match_y;temp];
    temp1 = unique(temp(:,1));
    for j = 1:length(temp1)
        if nnz(temp(:,1)==temp1(j)) > 1
            temp2 = temp(temp(:,1)==temp1(j),2);
            temp2 = unique(temp2);
            match_y =[match_y; [temp2(1:end-1),repmat(temp2(end),length(temp2)-1,1)]];
        end
    end
end
for i = 1:size(match_y,1)
    i
    Istitched(Istitched==match_y(i,2)) = match_y(i,1);
end

clear temp i j
csvwrite('Stitched\match_x_pass1.csv',match_x);
csvwrite('Stitched\match_y_pass2.csv',match_y);

clear match_x match_y

%% renumerate objects (optional, takes long time)
% [uni_cell,~,Irenumbered] = unique(Istitched(:));
% Irenumbered = uint32(reshape(Irenumbered-1,size(Istitched,1),[]));
uni_cell = unique(Istitched(:));
Irenumbered = Istitched;
clear Istitched
% save('Stitched\RenumberedCells.mat','Irenumbered_ExpCells','-v7.3','-append');

save('Stitched\RenumberedCells.mat','Irenumbered','-v7.3');
fid = fopen('Stitched\Cells_unique.txt','w');
fprintf(fid,'%d\n',uni_cell);
fclose(fid);

%% find "cells"
exp_d = 20;

[D,idx] = bwdist(Irenumbered);
D = uint32(D<=exp_d);

L = reshape(Irenumbered(idx(:)),size(Irenumbered));
CellLabel = D.*L;
% CellLabel = Irenumbered;
S.(['CellLabel_' num2str(exp_d)]) = CellLabel;
save('Stitched\RenumberedCells.mat','-struct','S','-append');

Outlines = CellLabel ~= imerode(CellLabel,strel('disk',2));
imwrite(Outlines,['Stitched\CellOutlines_' num2str(exp_d) '.png']);
% 
% back = imread('E:\PROOOJECTS\12_Neuron_mapping\160211_SamplePrep2\brain4_1\160211_UCL_b1_mip_100%_s7c1.jpg');
% % Outlines = imresize(Outlines,1);
% back = back + repmat(uint8(Outlines(1:size(back,1),1:size(back,2))*255),1,1,3);
% imwrite(back,'Stitched\DAPI_w_CellOutlines.png');

%% label child blobs
% blob = importdata('E:\PROOOJECTS\12_Neuron_mapping\160211_SamplePrep2\brain4_1\CP_160222\beforeQT_details.csv');
% blobspos = blob.data(:,1:2);
% blobspos_f = floor(blobspos+1);
% blobspos_idx = (blobspos_f(:,1)-1)*tile_size*tiles_y + blobspos_f(:,2);
% parentcell = CellLabel(blobspos_idx);
% csvwrite('ParentCell_beforeQT.csv',parentcell);

blob = importdata('..\CP_161014\Decoding\QT_0.45_1e-05_details.csv');
blobspos = blob.data(:,1:2);
blobspos_f = floor(blobspos+1);
blobspos_idx = (blobspos_f(:,1)-1)*tile_size*tiles_y + blobspos_f(:,2);
parentcell = CellLabel(blobspos_idx);
csvwrite('ParentCell\ParentCell_QT_0.45.csv',parentcell);

fid = fopen('ParentCell\QT_0.45_details_wCell_.csv','w');
formspec = [];
for i = 1:size(blob.textdata,2)+1
    formspec = [formspec, '%s,'];
end
formspec = [formspec(1:end-1), '\n'];
header = [blob.textdata(1,:),{'Parent_cell'}];
fprintf(fid,formspec,header{:});
blobwrite = [blob.textdata(2:end,1:2),num2cell([blob.data,double(parentcell)])]';
formspec = '%s,%s,';
for i = 1:size(blob.textdata,2)-1
    formspec = [formspec, '%d,'];
end
formspec = [formspec(1:end-1), '\n'];
fprintf(fid,formspec,blobwrite{:});
fclose(fid);

centroid = regionprops(CellLabel,'centroid');
centroid = cat(1,centroid.Centroid);
% csvwrite('Stitched\Cell_centroid.csv',centroid);
area = regionprops(CellLabel,'area');
area = cat(1,area.Area);
fid = fopen('Stitched\ExpCells.csv', 'w');
fprintf(fid, 'CellID,metadata_position,area,global_x_pos,global_y_pos\n');
fprintf(fid, '%d,%d,%d,%d,%d\n',...
    reshape([(1:length(area))',nan(length(area),1),area,centroid]',[],1));
fclose(fid);

