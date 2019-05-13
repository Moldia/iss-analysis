function [plotx_pool,ploty_pool,plotx_query,ploty_query,plotx_inter,ploty_inter,pool_id] = ...
    pairintersection_f(name_p1,name_p2,idx_re,pos,radius)

% Calculate the intersecting parts of two transcripts
% Xiaoyan 2014-11-18


%% position of the genes
pair1_x = pos(idx_re==name_p1,1);
pair1_y = pos(idx_re==name_p1,2);
pair2_x = pos(idx_re==name_p2,1);
pair2_y = pos(idx_re==name_p2,2);

%% re-order pairs
if length(pair1_x)>=length(pair2_x)
    pool_id = 1;
    pool = [pair1_x,pair1_y];
    query = [pair2_x,pair2_y];
    pool_x = pair1_x;
    pool_y = pair1_y;
    query_x = pair2_x;
    query_y = pair2_y;
else
    pool_id = 2;
    pool = [pair2_x,pair2_y];
    query = [pair1_x,pair1_y];
    pool_x = pair2_x;
    pool_y = pair2_y;
    query_x = pair1_x;
    query_y = pair1_y;
end

[plotx_pool,ploty_pool] = dot2poly_f(pool_x,pool_y,radius,15);
[plotx_query,ploty_query] = dot2poly_f(query_x,query_y,radius,15);

%% find intersections
[I,~] = rangesearch(pool,query,2*radius);   % first find neighbors within a fixed radius to save computational time

NN_num = zeros(length(I),1);
for i = 1:length(I)
    NN_num(i) = size(I{i},2);
end

idx_query = find(NN_num);
plotx_inter = {}; ploty_inter = {};

for i = 1:length(idx_query)
    idx_pool = I{idx_query(i)};
    for j = 1:length(idx_pool)
        [x,y] = polybool('intersection',...
            plotx_query(:,idx_query(i)),ploty_query(:,idx_query(i)),...
            plotx_pool(:,idx_pool(j)),ploty_pool(:,idx_pool(j)));
        
        if isempty(x)
        else
            plotx_inter = [plotx_inter; {x}];
            ploty_inter = [ploty_inter; {y}];
        end
    end
end

end
