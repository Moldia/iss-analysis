function I = resizestitch(ntiles, tilesize, resizef, immatrix)
% resize and stitch tile images
% Xiaoyan, 2017


% I = zeros(ntiles(2)*tilesize*resizef, ntiles(1)*tilesize*resizef, 3, 'uint8');
I = zeros(ntiles(2)*tilesize*resizef, ntiles(1)*tilesize, 'uint32');
for i = 1:ntiles(2)
    for j = 1:ntiles(1)
        temp = imresize(immatrix{i,j}(:,:,:),resizef);
        I(tilesize*resizef*(i-1)+1:tilesize*resizef*i,...
            tilesize*resizef*(j-1)+1:tilesize*resizef*j,:) = temp;
    end
end

end
