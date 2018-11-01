imageFixed = imread('');
imageMoving = imread('');

% control point selection
[pointsFixed,pointsMoving] = cpselect(imageMoving,imageFixed,'wait',true)