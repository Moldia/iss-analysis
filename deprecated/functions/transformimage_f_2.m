function original = transformimage_f_2(original,ref,yup,xleft,name,new_pos,scale)

if yup>=0
    original = original(yup+1:end,:);
else
    original = [zeros(-yup,size(original,2));original];
end
if xleft>=0
    original = original(:,xleft+1:end);
else
    original = [zeros(size(original,1),-xleft),original];
end

% remove extra columns/rows in the right lower corner so that all images
% have the same size
refsize = imfinfo(ref);
refsize = [refsize.Height,refsize.Width];
if refsize(1)<size(original,1)
    original = original(1:refsize(1),:);
else
    original = [original;zeros(refsize(1)-size(original,1),size(original,2))];
end
if refsize(2)<size(original,2)
    original = original(:,1:refsize(2));
else
    original = [original,zeros(size(original,1),refsize(2)-size(original,2))];
end
imwrite(original,name,'tiff','compression','none');

close all;
imshow(original);
hold on;
plot(new_pos(:,1)*scale,new_pos(:,2)*scale,'.')