function [connected_posx,connected_posy] = connectpoly(L,poly_x,poly_y)
% merge polygons more efficiently
% Xiaoyan 2014-11-13

connected_posx = {};
connected_posy = {};

for i = 1:length(L)
    if isempty(L{i})
    else
        connected_dot = L{i};
%         connected_dot = sort(connected_dot);
        connected_x = [];
        connected_y = [];
        
        % divide into groups with 10 transcripts
        walks = ceil(length(connected_dot)/10);
        for j = 1:walks
            if j == walks
                endidx = length(connected_dot);
            else
                endidx = 10*j;
            end
            temp_x = [];
            temp_y = [];
            for k = 10*(j-1)+1:endidx
                temp_x_back = temp_x;
                temp_y_back = temp_y;
                
                % polygon boolean operation - AND
                [temp_x,temp_y] = polybool('union',...
                    temp_x,temp_y,...
                    poly_x(:,connected_dot(k)),poly_y(:,connected_dot(k)));
                
                if isnan(polyarea(temp_x,temp_y))
                    temp_x = temp_x(~isnan(temp_x));
                    temp_y = temp_y(~isnan(temp_y));
                    if isnan(polyarea(temp_x,temp_y))
                        return
                    end
                end
                
            end

            % combine progressively
            if walks == 1
                connected_x = temp_x;
                connected_y = temp_y;
            else
%                 j*10
                [temp_x,temp_y] = poly2cw(temp_x,temp_y);
%                 [connected_x,connected_y] = poly2cw(connected_x,connected_y);
                [connected_x,connected_y] = polybool('union',...
                    temp_x,temp_y,connected_x,connected_y);
                connected_x = connected_x(~isnan(connected_x));
                connected_y = connected_y(~isnan(connected_y));
            end
            
        end
        
        connected_posx = [connected_posx;{connected_x}];
        connected_posy = [connected_posy;{connected_y}];
    end
end

end