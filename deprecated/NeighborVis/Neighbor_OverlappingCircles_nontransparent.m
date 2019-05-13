%% Overlapping circles to visualize neighbors
%  pair-wise, non-tranparent circles
%  Xiaoyan 2015-1-20

%%
format compact
warning('off','all');

%% parameters
decoding_file = 'input_example\DECODE_0.4_0.007_beforeQT_details_newform.csv';
pair = {'COL3A1' 'HER2'};
radius = 300;

%% transcripts
[name,pos] = getinsitudata_f(decoding_file);
 
% unique transcripts
[name_uni,~,idx_re] = unique(name);

% pair gene index
name_p1 = find(strcmp(name_uni,pair{1}));
name_p2 = find(strcmp(name_uni,pair{2}));
if isempty(name_p1) || isempty(name_p2)
    error('At least one of the genes specified does not have any positional information.');
end

%% images and plot
tic
x1 = pos(idx_re==name_p1,1)/5;
y1 = pos(idx_re==name_p1,2)/5;
x2 = pos(idx_re==name_p2,1)/5;
y2 = pos(idx_re==name_p2,2)/5;
Itemp1 = accumarray(ceil([y1,x1]),1,ceil([max(pos(:,2)),max(pos(:,1))]/5)+[1,1]);
Itemp2 = accumarray(ceil([y2,x2]),1,ceil([max(pos(:,2)),max(pos(:,1))]/5)+[1,1]);
Itemp1 = imdilate(Itemp1,strel('disk',radius/5));
Itemp2 = imdilate(Itemp2,strel('disk',radius/5));
I = cat(3,Itemp1,Itemp2,zeros(size(Itemp1)));
figure,imshow(I)
toc

title([pair{1} ' - ' pair{2} ' (r=' num2str(radius) ')']);


