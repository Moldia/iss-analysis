function new_pos = transform_positions_f(pos,scale,Isize,ang,yup,xleft,flo_rotate_size,fact)

%  Xiaoyan, 2015-6-26

pos_re = pos*scale;
pos_transformed = bsxfun(@minus,pos_re(:,1:2),floor(fliplr(Isize*scale)/2));

% rotation
rot_angle = -1*ang/180*pi;
rot_mat=[cos(rot_angle),sin(rot_angle);...
    -sin(rot_angle),cos(rot_angle)];
new_pos = pos_transformed*rot_mat;
rotate_size = flo_rotate_size*fact;

new_pos = bsxfun(@plus,fliplr(rotate_size)/2,new_pos);

% translation
new_pos(:,1) = new_pos(:,1)-xleft*fact;
new_pos(:,2) = new_pos(:,2)-yup*fact;
