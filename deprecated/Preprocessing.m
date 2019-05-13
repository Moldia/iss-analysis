%  tophat, then MIP

mkdir('MIP_tophat');

emptyt = csvread('EmptyTiles.csv', 1);
grid_x = 16;
grid_y = 12;
mmax = 168;
emptyt = (emptyt(:,2)-1)*grid_x + emptyt(:,1);
Iempty = zeros(2048,2048,'uint16');

for b = 1:3
    if b == 2 
        folder = ['161101_4028_36-2_b' num2str(b)];
        for c = 1:6
            if c == 1;
                SE = strel('disk', 10);
            else
                SE = strel('disk', 2);
            end
            
            t = 1;
            for m = 1:mmax
                % originally empty FOV
                while ismember(t, emptyt)
                    tnum = num2str(t);
                    while length(tnum) < 3
                        tnum = ['0' tnum];
                    end
                    imwrite(Iempty, ['MIP_tophat\base' num2str(b) '_c' num2str(c) 't' tnum '_MIP.tif'], 'tiff');
                    t = t+1;
                end
                    
                mnum = num2str(m);
                while length(mnum) < 3
                    mnum = ['0' mnum];
                end

%                 iminfo = imfinfo([folder '\' folder '_z1c1m' mnum '_ORG.tif']);
                MIP = zeros(2048, 2048, 'uint16');

                for z = 1:7
                    fname = [folder '\' folder '_z' num2str(z) 'c' num2str(c) 'm' mnum '_ORG.tif'];
                    I = imread(fname);

                    % tophat
                    Itop = imtophat(I, SE);
                    MIP = max(MIP, Itop);

                end

                tnum = num2str(t);
                while length(tnum) < 3
                    tnum = ['0' tnum];
                end
                imwrite(MIP, ['MIP_tophat\base' num2str(b) '_c' num2str(c) 't' tnum '_MIP.tif'], 'tiff');
                t = t+1;
            end
            
            while t <= grid_x*grid_y
                tnum = num2str(t);
                while length(tnum) < 3
                    tnum = ['0' tnum];
                end
                imwrite(Iempty, ['MIP_tophat\base' num2str(b) '_c' num2str(c) 't' tnum '_MIP.tif'], 'tiff');
                t = t+1;
            end
        end
    end
end


    