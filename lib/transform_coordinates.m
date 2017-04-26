function new_pos = transform_coordinates(pos, scale, imsize, ang, yup, xleft,...
    imsize_final, imscale)
% Xiaoyan, 2017

% scaling first
pos_re = pos*scale;

% center as origo
pos_transformed = bsxfun(@minus, pos_re(:,1:2), floor(fliplr(imsize*scale)/2));

% rotation
rot_angle = -1*ang/180*pi;
rot_mat = [...
    cos(rot_angle), sin(rot_angle);...
    -sin(rot_angle), cos(rot_angle)];
new_pos = pos_transformed*rot_mat;

posscale = 1/imscale;
rotate_size = imsize_final*posscale;

new_pos = bsxfun(@plus, fliplr(rotate_size)/2, new_pos);

% translation
new_pos(:,1) = new_pos(:,1)-xleft*posscale;
new_pos(:,2) = new_pos(:,2)-yup*posscale;

end
