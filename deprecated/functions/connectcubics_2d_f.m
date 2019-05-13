function I = connectcubics_2d_f(grid_size,Isize,pos,up,down,zposition,fh,col)
% connect neighboring cubics for better visualization
% Xiaoyan, 2014-12-2

grid_nrx =  ceil(Isize(2)/grid_size);
grid_nry =  ceil(Isize(1)/grid_size);

transformed_x = ceil(pos(:,1)/grid_size);
transformed_y = ceil(pos(:,2)/grid_size);

I = accumarray([transformed_y,transformed_x],1,[grid_nry,grid_nrx]);

if isnumeric(fh)
    fh = figure; hold on;
end
for j = 1:grid_nrx
    for i = 1:grid_nry
        if I(i,j)
            wall_nr = zeros(6,1);
            if i == 1
                wall_nr(1) = 1;
                if j == 1
                    wall_nr(2) = 1;
                elseif j == grid_nrx
                    wall_nr(3) = 1;
                else
                    if ~I(i,j+1)
                        wall_nr(4) = 1;
                    end
                    if ~I(i,j-1)
                        wall_nr(2) = 1;
                    end
                end
            else
                if j == 1
                    wall_nr(2) = 1;
                elseif j == grid_nrx
                    wall_nr(3) = 1;
                else
                    if ~I(i-1,j)
                        wall_nr(1) = 1;
                    end
                    if ~I(i,j-1)
                        wall_nr(2) = 1;
                    end
                    if ~I(i+1,j)
                        wall_nr(3) = 1;
                    end
                    if ~I(i,j+1)
                        wall_nr(4) = 1;
                    end
                end
            end
            if up
                wall_nr(5) = 1;
            end
            if down
                wall_nr(6) = 1;
            end
            drawwalls_f(grid_size,...
                [j-1 i-1 zposition].*[grid_size grid_size 1],...
                wall_nr,col,0.2,'none',fh);
        end
    end
end
                        
            
        