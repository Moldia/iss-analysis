%% A and B color distribution
%  simulations
%  Xiaoyan, 2015-6-5

%% parameters - simulation
imsize = 200;
num_blobs_1 = 500;
ratio_blobs = 0.7;  % between 0 and 1
distance = 2;

%% blob 1 pseudorandom distribution
randx = randi(imsize,num_blobs_1,1);
randy = randi(imsize,num_blobs_1,1);
blobs_1 = accumarray([randy,randx],1,[imsize,imsize]);
c = 0;
while nnz(blobs_1>=2)
    c = c+1;
    if c >= 10
        error('too many blobs!')
    else
        nonzero = nnz(blobs_1>=2);
        blobs_1(blobs_1>=2) = 1;
        randx = randi(imsize,nonzero,1);
        randy = randi(imsize,nonzero,1);
        blobs_1 = blobs_1 + accumarray([randy,randx],1,[imsize,imsize]);
    end
end
label_1 = find(blobs_1);

subplot(2,2,1), imshow(blobs_1);
% [m,n] = find(blobs_1);
% Dist = pdist([m,n]);
dist_blobs_1 = bwdist(blobs_1);
subplot(2,2,2), imshow(max(dist_blobs_1(:))-dist_blobs_1,[]);

Water_1 = watershed(dist_blobs_1);
L = label2rgb(Water_1);
subplot(2,2,3), imshow(L)

%% blob 2 distribtuions
num_blobs_2 = round(num_blobs_1*ratio_blobs,0);
num_dist = nnz(dist_blobs_1==distance);
if max(Water_1(:))<num_blobs_2
    error('too high density');
else
%     blobs_2_idx = datasample(find(dist_blobs_1>=distance),num_blobs_2); % random sampling
%     blobs_2 = zeros(imsize^2,1);
%     blobs_2(blobs_2_idx) = 1;
%     blobs_2 = reshape(blobs_2,imsize,[]);
    blobs_2_1 = double(datasample(1:max(Water_1(:)),num_blobs_2,'replace',false));
    blobs_2_pos = zeros(length(blobs_2_1),1);
    Pool_extra = [];
    for i = 1:length(blobs_2_1)
        parent = find(Water_1==blobs_2_1(i) & dist_blobs_1==distance);
        pool_extra = [];
        while isempty(parent)
            pool_temp = 1:max(Water_1(:));
            pool_temp = ~ismember(Water_1,blobs_2_1) & ~ismember(Water_1,Pool_extra);
            if isempty(pool_temp)
                error('no more possible non-overlapping sampling');
            else
                pool_extra = datasample(Water_1(pool_temp),1);
                parent = find(Water_1==pool_extra & dist_blobs_1==distance);
            end
        end

        blobs_2_pos(i) = datasample(parent,1);
        Pool_extra = [Pool_extra,pool_extra];
    end
        
    blobs_2 = zeros(imsize^2,1);
    blobs_2(blobs_2_pos) = 1;
    blobs_2 = reshape(blobs_2,imsize,[]);        
end
% subplot(2,2,3), imshow(blobs_2);

subplot(2,2,4), imshow(cat(3,blobs_1,blobs_2,blobs_1));
    