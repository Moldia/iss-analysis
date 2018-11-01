
%% parameters
[names,Pos] = getinsitudata_f('E:\PROOOJECTS\11_PDGFR\K\CP_150708\Decoding\QT_0.4_0.001_details.csv');
[name_uni,~,idx_re] = unique(names);
background = imread('E:\PROOOJECTS\11_PDGFR\K\IHC\IHC_mergedRGB_20%.tif');
scale = 0.2;

name_density = 'CDH1';
max_dist = 700;

%% DBSCAN
idx_density = find(strcmp(name_uni,name_density));

if isempty(idx_density)
    error('No specified transcript detected in the input file');
end
pos_density = Pos(idx_re==idx_density,1:2);
[lab,labc] = dbscan(pos_density,max_dist,20);

pos_density = correctcoordinates_f(pos_density,scale);

figure,imshow(background),hold on;
for i = 1:max(lab)
    plot(pos_density(lab==i,1),pos_density(lab==i,2),'+');
end

