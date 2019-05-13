function new_pos = transformpositions_f(matfile,scale,ang,yup,xleft,flo_rotate_size,fact)

%  Xiaoyan, 2015-5-26


load(matfile);
pos_re(:,1:2) = pos_re(:,1:2)*scale;
pos_transformed = bsxfun(@minus,pos_re(:,1:2),floor(fliplr(Isize*scale)/2));

% rotation
rot_angle = -1*ang/180*pi;
rot_mat=[cos(rot_angle),sin(rot_angle);...
    -sin(rot_angle),cos(rot_angle)];
new_pos = pos_transformed*rot_mat;
% rotate_size = max([size(flo_rotate);size(ref)],[],1)*5;
rotate_size = flo_rotate_size*fact;

new_pos = bsxfun(@plus,fliplr(rotate_size)/2,new_pos);

% translation
new_pos(:,1) = new_pos(:,1)-xleft*fact;
new_pos(:,2) = new_pos(:,2)-yup*fact;

new_pos = [new_pos,pos_re(:,3:4)];
