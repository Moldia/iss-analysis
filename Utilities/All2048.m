
for b = 1:4
    mkdir(['All2048\base' num2str(b)]);
end

% make sure all images are exactly 2048x2048 in size
for b = 1:4
    for s = 1:28
        if b == 3
            sourceprefix = fullfile(['Base' num2str(b)],...
            ['new-Orthogonal Projection-01_s' paddigits(s, 2) 'c']);
        else
            sourceprefix = fullfile(['Base' num2str(b)],...
            ['IP-Orthogonal Projection-01_s' paddigits(s, 2) 'c']);
        end
        
		% get image file information
        imsize = imfinfo([sourceprefix '1_ORG.tif']);
        imsize = [imsize.Width, imsize.Height];
        
        if isequal(imsize, [2048 2048])
            for c = 1:6
                source = [sourceprefix num2str(c) '_ORG.tif'];
                destination = fullfile('All2048', ['base' num2str(b)],...
                    ['base' num2str(b) '_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif']);
                copyfile(source, destination);
            end
        else
            for c = 1:6
                source = [sourceprefix num2str(c) '_ORG.tif'];
                destination = fullfile('All2048', ['base' num2str(b)],...
                    ['base' num2str(b) '_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif']);
                
                I = imread(source);
                I = padimg(I, 2048-imsize(1), 2048-imsize(2));
                imwrite(I, destination);
            end
        end
    end
end

