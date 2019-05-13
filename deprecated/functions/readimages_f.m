function [ref,flo,original,Ifuse] = readimages_f(ref,isrgbref,scaleref,flo,isrgbflo,scaleflo)

ref = imread(ref);
flo = imread(flo);

if isrgbref
    ref = ref(:,:,3);
end
if isrgbflo
    flo = flo(:,:,3);
end

original = imresize(flo,0.4);

ref = imresize(ref,scaleref);
flo = imresize(flo,scaleflo);

Ifuse = imfuse(flo,ref);
