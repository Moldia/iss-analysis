%% visualization
%  data source: E:\PROOOJECTS\9_Lung_Multiplex\8128_4\CPresult_140926_NucleiAlign_ExpandedBlobs
%  22/10/2014 Xiaoyan


%% load decoded result
load('E:\PROOOJECTS\9_Lung_Multiplex\8128\8128_4\CPresult_140926_NucleiAlign_ExpandedBlobs\Decoding\blob.mat')

measurements = csvread('E:\PROOOJECTS\9_Lung_Multiplex\8128\8128_4\CPresult_140926_NucleiAlign_ExpandedBlobs\ExpandedBlobs.csv',1);
pos = measurements(:,14:15);

clear mearuements

%% create an image representing the tissue
I = imread('E:\PROOOJECTS\9_Lung_Multiplex\8128\8128_4\ZENout\8128_4_c1_ORG_padded.tif');
I_modified = im2bw(I,graythresh(I));
I_modified = imdilate(I_modified,strel('disk',50));
I_modified = imfill(I_modified,'holes');
% imshow(I_modified)

%% trace boundaries
% boundaries = bwboundaries(I_modified);
% hold on;
% for k=1:10
%    b = boundaries{k};
%    plot(b(:,2),b(:,1),'g','LineWidth',3);
% end

%% divide into grids
grid_size = 400;

grid_nrx =  ceil(size(I,2)/grid_size);
grid_nry =  ceil(size(I,1)/grid_size);
grid_nr = grid_nrx*grid_nry;

%% background grid and 3-d representation
I_modified = [I_modified,zeros(size(I_modified,1),grid_nrx*grid_size-size(I_modified,2))];
I_modified = [I_modified;zeros(grid_nry*grid_size-size(I_modified,1),size(I_modified,2))];

polygoncoord = zeros(5,2,grid_nr);
grid_background = zeros(grid_nry,grid_nrx);
% figure; hold on
for j = 1:grid_nrx   % j along x axis, column
    gridxmin = (j-1)*grid_size;
    gridxmax = j*grid_size;
    
    for i = 1:grid_nry  % i along y axis, row
        gridymin = (i-1)*grid_size;
        gridymax = i*grid_size;
        
        k = (j-1)*grid_nry+i; % counting direction: y
        
        polyx = [gridxmin,gridxmin,gridxmax,gridxmax,gridxmin];
        polyy = [gridymin,gridymax,gridymax,gridymin,gridymin];
        polygoncoord(:,:,k) = [polyx',polyy'];
        
        if length(find(I_modified(gridymin+1:gridymax-1,gridxmin+1:gridxmax-1)))>= .3*grid_size^2
            grid_background(i,j) = 1;
            center = [(gridxmin+gridxmax)/2,(gridymin+gridymax)/2,0];
            rppd([340 340 340],center,rgb('lightskyblue'),[1 1 1],.05);
            center = [(gridxmin+gridxmax)/2,(gridymin+gridymax)/2,200];
            rppd([340 340 340],center,rgb('lightskyblue'),[1 1 1],.05);
            center = [(gridxmin+gridxmax)/2,(gridymin+gridymax)/2,600];
            rppd([340 340 340],center,rgb('lightskyblue'),[1 1 1],.05);

        end
        
    end
end
clear gridxmin gridxmax gridymin gridymax polyx polyy
set(gca,'YDir','reverse',...
    'PlotBoxAspectRatio',[1 size(I_modified,1)/size(I_modified,2) size(I_modified,1)/1500]);
% axis equal
axis off

%% visualize grid size
% imshow(I_modified)
% hold on;
% plot(polygoncoord(:,1,floor(grid_nrx/2)*grid_nry+floor(grid_nry/2)),...
%     polygoncoord(:,2,floor(grid_nrx/2)*grid_nry+floor(grid_nry/2)),'b');

%% ACTB test
% figure,hold on
ACTB_in = zeros(grid_nr,1);
temp_in = [];
for j = 1:grid_nrx   % j along x axis, column    
    for i = 1:grid_nry  % i along y axis, row
        
        k = (j-1)*grid_nry+i; % counting direction: y
        
        polyx = polygoncoord(:,1,k);
        polyy = polygoncoord(:,2,k);
        
        temp_in =  inpolygon(pos(blob==1,1),pos(blob==1,2),polyx,polyy);
        
        ACTB_in(k) = length(find(temp_in));
       
    end
end

ACTB_in = ACTB_in/max(ACTB_in);

for j = 1:grid_nrx   % j along x axis, column
    gridxmin = (j-1)*grid_size;
    gridxmax = j*grid_size;
    
    for i = 1:grid_nry  % i along y axis, row
        gridymin = (i-1)*grid_size;
        gridymax = i*grid_size;
        
        k = (j-1)*grid_nry+i;
        
        center = [(gridxmin+gridxmax)/2,(gridymin+gridymax)/2,1+grid_size*2];
        if ACTB_in(k) == 0
            alph = 0;
        else
            alph = ceil(ACTB_in(k)*10)/10;
        end
        
        if grid_background(i,j) == 1;
            rppd([340 340 340],center,rgb('black'),[1 1 1],alph);
        end

        
    end
end

set(gca,'YDir','reverse',...
    'PlotBoxAspectRatio',[1 size(I_modified,1)/size(I_modified,2) 20/340]);


%% AF750, Cy3, Cy5 signals
signals_in = zeros(grid_nr,3);
for j = 1:grid_nrx   % j along x axis, column    
    for i = 1:grid_nry  % i along y axis, row
        
        k = (j-1)*grid_nry+i; % counting direction: y
        
        polyx = polygoncoord(:,1,k);
        polyy = polygoncoord(:,2,k);

        for c = 1:3
            temp_in =  inpolygon(pos(blob==c,1),pos(blob==c,2),polyx,polyy);
            signals_in(k,c) = length(find(temp_in));
        end
       
    end
end

for c = 1:3
    signals_in(:,c) = signals_in(:,c)/max(signals_in(:,c));
end
signals_in(:,[1,3]) = zeros(size(signals_in,1),2);

for j = 1:grid_nrx   % j along x axis, column
    gridxmin = (j-1)*grid_size;
    gridxmax = j*grid_size;
    
    for i = 1:grid_nry  % i along y axis, row
        gridymin = (i-1)*grid_size;
        gridymax = i*grid_size;
        
        k = (j-1)*grid_nry+i;
        
%         center = [(gridxmin+gridxmax)/2,(gridymin+gridymax)/2,1+grid_size*4];
        center = [(gridxmin+gridxmax)/2,(gridymin+gridymax)/2,200];
%         if ACTB_in(k) == 0
%             alph = 0;
%         else
%             alph = ceil(ACTB_in(k)*10)/10;
%         end
%         
        if grid_background(i,j) == 1;
            if signals_in(k,:) == [0 0 0]
                rppd([340 340 340],center,[1 1 1],[1 1 1],.2);
            else
                rppd([340 340 340],center,ceil(signals_in(k,:)),[1 1 1],.2);
            end
        end

        
    end
end

% set(gca,'YDir','reverse',...
%     'PlotBoxAspectRatio',[1 size(I_modified,1)/size(I_modified,2) 20*5/340]);
% axis off
