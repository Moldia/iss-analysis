% image registration using DFT
% Xiaoyan, 2017

%% DFT registration first

for b = 1:4
    mkdir(['base' num2str(b)]);
end

% base 1 make sure all images are exactly 2048x2048 in size
for s = 1:28
    imsize = imfinfo(['..\Base1\new-Orthogonal Projection-01_s' paddigits(s, 2) 'c1_ORG.tif']);
    imsize = [imsize.Width, imsize.Height];
    if isequal(imsize, [2048 248])
        for c = 1:6
            source = ['..\Base1\new-Orthogonal Projection-01_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif'];
            destination = ['base1\base1_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif'];
            copyfile source destination
        end
    else
        for c = 1:6
            source = ['..\Base1\new-Orthogonal Projection-01_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif'];
            destination = ['base1\base1_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif'];
            I = imread(source);
            I = padimg(I, 2048-imsize(1), 2048-imsize(2));
            imwrite(I, destination);
        end
    end
end


%% the rest align to base1
fails = [];
for s = 1:28
    for b = 2:4
        % registration on general stain
        ref = imread(['base1\base1_s' paddigits(s, 2) 'c2_ORG.tif']);
        float = imread(fullfile(['..\Base' num2str(b)],...
            ['IP-Orthogonal Projection-01_s' paddigits(s, 2) 'c2_ORG.tif']));
        
        % pad if image sizes don't match
        if ~isequal(size(float), [2048 2048])
            float = padimg(float, 2048-size(float,2), 2048-size(float,1));
        end
        
        % registration resolution 0.1 px
        [tform, Greg] = dftregistration(fft2(ref), fft2(float), 10);
        
        % use results only if error is < 0.12
        if tform(1) < 0.12
            
            % inverse FFT of floating image
            float2 = ifft2(Greg);
            float2 = uint16(real(float2));
            
            % write
            destination = fullfile(['base' num2str(b)],...
                ['base' num2str(b) '_s' paddigits(s, 2) 'c2_ORG.tif']);
            imwrite(float2, destination);
            
            % process rest of the channels
            for c = 1:6
                if c ~= 2
                    float = imread(['..\Base2\IP-Orthogonal Projection-01_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif']);
                    if ~isequal(size(float), [2048 2048])
                        float = padimg(float, 2048-size(float,2), 2048-size(float,1));
                    end
                    
                    % use the same tform
                    [nr,nc] = size(float);
                    Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
                    Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
                    
                    [Nc,Nr] = meshgrid(Nc,Nr);
                    Greg = fft2(float).*exp(1i*2*pi*(-tform(3)*Nr/nr-tform(4)*Nc/nc));
                    Greg = Greg*exp(1i*tform(2));
                    
                    % inverse FFT
                    float2 = ifft2(Greg);
                    float2 = uint16(real(float2));
                    
                    % write
                    destination = fullfile(['base' num2str(b)],...
                        ['base' num2str(b) '_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif']);
                    imwrite(float2, destination);
                end
            end
        else
            fails = [fails; s, b];
            for c = 1:6
                source = ['..\Base2\IP-Orthogonal Projection-01_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif'];
                destination = fullfile(['base' num2str(b)],...
                        ['base' num2str(b) '_s' paddigits(s, 2) 'c' num2str(c) '_ORG.tif']);
                copyfile(source, destination);
            end
        end
    end
end

figure;
for i = 1:12
    subplot(3,4,i);
    I1 = imread(['base1\base1_s' paddigits(num2str(fails(i,1)),2) 'c2_ORG.tif']);
    I2 = imread(['base' num2str(fails(i,2)) '\base' num2str(fails(i,2)) '_s' paddigits(num2str(fails(i,1)),2) 'c2_ORG.tif']);
        
    imshowpair(I1, I2);
    title(['base' num2str(fails(i,2)) ' scene' num2str(fails(i,1))]);
end

        
         