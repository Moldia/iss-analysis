function Isections = combinesections_f(data,grid_size,Isize)
% combine data from consecutive sections and make it into 3D matrix
% Xiaoyan, 2014


grid_size = 50;

grid_nrx =  ceil(Isize(2)/grid_size);
grid_nry =  ceil(Isize(1)/grid_size);

Isections = zeros(grid_nry,grid_nrx,length(data));

for i = 1:length(data)
    pos = data{i};
    transformed_x = ceil(pos(:,1)/grid_size);
    transformed_y = ceil(pos(:,2)/grid_size);
    Itemp = accumarray([transformed_y,transformed_x],1,[grid_nry,grid_nrx]);
    Itemp = Itemp/sum(Itemp(:))*500;

    Isections(:,:,i) = Itemp;
end


