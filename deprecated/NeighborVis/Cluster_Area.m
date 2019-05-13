%% Neighbor area
%  pair-wise
%  also see Neighbor_OverlappingCircles

%%
format compact
warning('off','all');

%% parameters
seq = importdata('input_example\QT_0.3_0.01_details.csv',',',1);
pair = {'PDGFRA' 'PDGFRB'};
image = imread('input_example\RGB_adjusted3.jpg');
radius = 200;

%% extract data
gene_name = seq.textdata(2:end,2);
position = seq.data(:,1:2);
clear seq;
 
% unique transcripts
[uniname,~,uniname_re] = unique(gene_name);

% pair gene index
idx_p1 = find(strcmp(uniname,pair{1}));
idx_p2 = find(strcmp(uniname,pair{2}));
if isempty(idx_p1) || isempty(idx_p2)
    error('At least one of the genes specified does not have any positional information')
end

% position of the genes
pair1_x = position(uniname_re==idx_p1,1);
pair1_y = position(uniname_re==idx_p1,2);
pair2_x = position(uniname_re==idx_p2,1);
pair2_y = position(uniname_re==idx_p2,2);

%% re-order of pool and query
if length(pair1_x)>=length(pair2_x)
    pool_id = 1;
    pool_idx = idx_p1;
    query_idx = idx_p2;
    pool = [pair1_x,pair1_y];
    query = [pair2_x,pair2_y];
    pool_x = pair1_x;
    pool_y = pair1_y;
    query_x = pair2_x;
    query_y = pair2_y;
else
    pool_id = 2;
    pool_idx = idx_p2;
    query_idx = idx_p1;
    pool = [pair2_x,pair2_y];
    query = [pair1_x,pair1_y];
    pool_x = pair2_x;
    pool_y = pair2_y;
    query_x = pair1_x;
    query_y = pair1_y;
end

%% circles
alph = linspace(0,2*pi,21);
circx = radius*cos(alph);
circy = radius*sin(alph);
Area_solo = polyarea(circx,circy);

%% polygon coordinates
plotx_pool = (repmat(pool_x,1,21)+repmat(circx,length(pool_x),1))';
ploty_pool = (repmat(pool_y,1,21)+repmat(circy,length(pool_y),1))';
plotx_query = (repmat(query_x,1,21)+repmat(circx,length(query_x),1))';
ploty_query = (repmat(query_y,1,21)+repmat(circy,length(query_y),1))';

%% connect and merge polygons
[L1,NN_num] = List_self_f(position,query_idx,uniname_re,radius);
[connected_posx_query,connected_posy_query] = connectpoly_f(L1,plotx_query,ploty_query);
disp('conecting finished')

[L2,NN_num_2] = List_self_f(position,pool_idx,uniname_re,radius);
[connected_posx_pool,connected_posy_pool] = connectpoly_f(L2,plotx_pool,ploty_pool);
disp('conecting finished')

%% plot and area
figure; hold on;
imagesc(image); axis image;
axis off;
set(gca,'YDir','reverse');

Area_query = [];
for i = 1:length(connected_posx_query)
    patch(connected_posx_query{i},connected_posy_query{i},'r','edgecolor','b');
    Area_query = [Area_query;polyarea(connected_posx_query{i},connected_posy_query{i})];
end

solo = find(NN_num==0);
patch(plotx_query(:,solo),ploty_query(:,solo),'k','edgecolor','none');
Area_query = [Area_query;length(solo)*Area_solo];

% patch(plotx_query,ploty_query,'y','facealpha',.5,'edgecolor','none');

sum(Area_query,1)                

Area_pool = [];
for i = 1:length(connected_posx_pool)
    patch(connected_posx_pool{i},connected_posy_pool{i},'b','edgecolor','b');
    Area_pool = [Area_pool;polyarea(connected_posy_pool{i},connected_posy_pool{i})];
end

solo = find(NN_num_2==0);
patch(plotx_pool(:,solo),ploty_pool(:,solo),'c','edgecolor','none');
Area_query = [Area_query;length(solo)*Area_solo];

% patch(plotx_query,ploty_query,'y','facealpha',.5,'edgecolor','none');

sum(Area_query,1)  
                