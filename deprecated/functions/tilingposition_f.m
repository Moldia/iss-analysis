function tilingposition_f(img, scale, x_tile_size_original, y_tile_size_original)
% function to tiling position
% Xiaoyan, 2016-9-12

disp('plotting to show tiled position..');
Iref = imread(img);
imgsize = size(Iref);
Iref = imresize(Iref,0.1/scale);

tx = ceil(imgsize(2)/scale/x_tile_size_original);
ty = ceil(imgsize(1)/scale/y_tile_size_original);

plot_x = zeros(ty+1,tx+1);
plot_y = zeros(ty+1,tx+1);

for i=1:ty
    for j=1:tx
        tilenum = (i-1)*tx+j;
        plot_x(i+1,j+1) = x_tile_size_original*j;
        plot_y(i+1,j+1) = y_tile_size_original*i;
    end
end
plot_x(1,2:end) = plot_x(2,2:end);
plot_y(2:end,1) = plot_y(2:end,2);

figure; imagesc(Iref); axis image; hold on;
plot_x = floor(0.1*plot_x);
plot_y = floor(0.1*plot_y);

% vertical lines
%--------------------
for i = 1:tx+1
    x = plot_x(:,i);
    y = plot_y(:,i);
    plot (x,y,'y');
end
% horizental lines
%--------------------
for i = 1:ty+1
    x = plot_x(i,:);
    y = plot_y(i,:);
    plot (x,y,'y');
end
axis off;
end
