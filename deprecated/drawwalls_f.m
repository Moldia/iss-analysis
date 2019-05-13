function drawwalls_f(grid_size,start_point,wall_logic,col,alph,linestyl,fh)
% draw every wall severately for a cubic
% Xiaoyan, 2014-12-2

patch_back = [0,0,0;grid_size,0,0;grid_size,0,grid_size;0,0,grid_size;0,0,0];
patch_left = [0,0,0;0,grid_size,0;0,grid_size,grid_size;0,0,grid_size;0,0,0];
patch_front = [0,grid_size,0;grid_size,grid_size,0;grid_size,grid_size,grid_size;0,grid_size,grid_size;0,grid_size,0];
patch_right = [grid_size,grid_size,0;grid_size,0,0;grid_size,0,grid_size;grid_size,grid_size,grid_size;grid_size,grid_size,0];
patch_up = [0,0,grid_size;grid_size,0,grid_size;grid_size,grid_size,grid_size;0,grid_size,grid_size;0,0,grid_size];
patch_bottom = [0,0,0;grid_size,0,0;grid_size,grid_size,0;0,grid_size,0;0,0,0];

wall = [{patch_back},{patch_left},{patch_front},{patch_right},{patch_up},{patch_bottom}];


for i = 1:6
    if wall_logic(i)
        polygonwall = wall{i} + repmat(start_point,5,1);
        figure(fh);
        patch(polygonwall(:,1),polygonwall(:,2),polygonwall(:,3),col,...
            'facealpha',alph,'linestyle',linestyl);
    end
end