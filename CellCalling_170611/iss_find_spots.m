function [GoodGlobalYX, GoodSpotColors, GoodIsolated] = iss_find_spots(o)
% [GlobalYX, SpotColors, Isolated] = iss_find_spots(o)
%
% finds spots in all tiles using the reference channel, removes
% duplicates in overlap regions and returns nSpots x 2 array GlobalYX of
% coordinates in global frame
% 
% Looks up colors from apporpriate files and makes nSpots x nBP x nRounds
% array SpotColors
%
% final output Isolated is a nSpots x 1 binary array giving 1 for
% well-isolated spots
%
% NB spots that can't be read in all rounds are discarded
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
%% variable naming conventions:
% spot subgroups:
% All: Any spot included in any tile (includes duplicates)
% nd: only spots whose anchor coordinate is in its home tile (no duplicates)
% Good: Spots for which could be read for all rounds

% coordinate frames or other info
% LocalYX: relative to home tile on the reference round
% LocalTile: number of home tile on the reference round
% GlobalYX: relative to the stitched image on the reference round
% RoundYX: relative to home tile after registration on each round
% RoundTile: number of home tile after registration on each round
% Isolated: binary number, saying if it is isolated
% SpotColors: the answer:

%% basic variables
rr = o.ReferenceRound;
EmptyTiles = strcmp('', squeeze(o.TileFiles(rr,:,:)));
Tiles = find(~EmptyTiles)';

[nY, nX] = size(EmptyTiles);
nTiles = nY*nX;

%% loop through all tiles, finding spots in anchor channel on ref round
RawLocalYX = cell(nTiles,1);  % cell array, giving spots in local coordinates
RawIsolated = cell(nTiles,1);
for t=Tiles
    if mod(t,10)==0; fprintf('Detect spots in tile %d\n', t); end;
    [y,x] = ind2sub([nY nX], t);
    AnchorIm = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);
    [RawLocalYX{t}, RawIsolated{t}] = iss_detect_spots(AnchorIm, o);
end
    
%% now make array of global coordinates
AllIsolated = logical(vertcat(RawIsolated{:})); % I HATE MATLAB - for converting logical to doubles for no reason
nAll = length(AllIsolated);

AllGlobalYX = zeros(nAll,2);
AllLocalYX = zeros(nAll,2);

ind = 1;
for t=Tiles
    MySpots = RawLocalYX{t};
    nMySpots = size(MySpots, 1);
    AllGlobalYX(ind:ind+nMySpots-1,:) = bsxfun(@plus, MySpots, o.RefPos(t,:));
    AllLocalYX(ind:ind+nMySpots-1,:) = MySpots;
    OriginalTile(ind:ind+nMySpots-1) = t;
    ind = ind+nMySpots;
end
if o.Graphics
    figure(1001)
    plot(AllGlobalYX(:,2), AllGlobalYX(:,1), '.', 'markersize', 1);
    title('All global coords including duplicates');
    %set(gca, 'YDir', 'reverse');
end

%% now remove duplicates by keeping only spots detected on their home tile

[AllLocalTile, ~] = iss_which_tile(AllGlobalYX, o.RefPos, o.TileSz);
NotDuplicate = (AllLocalTile==OriginalTile');
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndLocalYX = AllLocalYX(NotDuplicate,:);
ndIsolated = AllIsolated(NotDuplicate,:);
ndLocalTile = AllLocalTile(NotDuplicate,:);

nnd = sum(NotDuplicate);

if o.Graphics
    figure(1002)
    plot(ndGlobalYX(:,2), ndGlobalYX(:,1), '.', 'markersize', 1);
    title('Global coords without duplicates');
    %set(gca, 'YDir', 'reverse');
end


%% find coordinates of each spot in appropriate tile
ndRoundYX = nan(nnd,2,o.nRounds);
ndRoundTile = nan(nnd,o.nRounds);

for r=1:o.nRounds
    fprintf('Finding coordinates for round %d\n', r);
    % for each spot, find which tile it is in for this round
    
    % this array contains the offset a single tile in round r to ref round
    % the last thing is a way of diagonalizing
    SameTileRelativePos = reshape(o.RelativePos(r,:,1:nTiles+1:nTiles^2), [2, nTiles])';
    
    [ndRoundTile(:,r), ~] = iss_which_tile(ndGlobalYX, o.RefPos-SameTileRelativePos, o.TileSz);
    
    % now shift it. sub2ind avoids making a 3d matrix. Maybe use IndexArrayNan?
    %EachSpotShift = squeeze(o.RelativePos(r,:,sub2ind([nTiles nTiles], ndRoundTile(:,r), ndLocalTile)))';
    if 0; %any(~cellfun(@isempty,o.PcTransform(r,:)))
        for t2=Tiles % tile you are coming from (home on ref round)
            for t1=Tiles % tile you are going to (on other round)
                MySpots = (ndRoundTile(:,r)==t1 & ndLocalTile==t2);
                if ~any(MySpots); continue; end;
                X1 = [ndLocalYX(MySpots,:), ones(sum(MySpots),1)];
                Y = clip(round(X1*o.PcTransform{r,t1,t2}),0,o.TileSz-1);
                ndRoundYX(MySpots,:,r) = Y;
            end
        end
    else
        IndexArray = zeros(4, nnd, 2);
        IndexArray(:,:,1) = [repmat([r 1], nnd, 1), ndRoundTile(:,r), ndLocalTile]';
        IndexArray(:,:,2) = [repmat([r 2], nnd, 1), ndRoundTile(:,r), ndLocalTile]';
        EachSpotShift = IndexArrayNan(o.RelativePos, IndexArray);
        ndRoundYX(:,:,r) = ndLocalYX + EachSpotShift;
    end
    
end

%% get spot colors
ndSpotColors = nan(nnd, o.nBP, o.nRounds);
ndCorrectedYX = nan(nnd,2,o.nRounds, o.nBP);

% MyTile = 42;
% MyNeighbors = MyTile + bsxfun(@plus, [-12 0 12], [-1 0 1]');
% for t=MyNeighbors(:)' 
for t=1:nTiles
    if EmptyTiles(t); continue; end

    if mod(t,10)==0; fprintf('reading spot colors for tile %d\n', t); end
    for r=1:o.nRounds
        MySpots = (ndRoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        for b=1:o.nBP
%           BaseIm = imread(FileName, o.AnchorChannel + b);
            TifObj.setDirectory(o.AnchorChannel + b);
            BaseIm = TifObj.read();
            if o.SmoothSize
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
            else
                BaseImSm = BaseIm;
            end
            if o.PointCloud
                % find spots for base b on tile t
                BaseYX = iss_detect_spots(BaseIm,o);
                % now loop over all potential home tiles for this one
                MyRefTiles = unique(ndLocalTile(ndRoundTile(:,r)==t));
                for t2 = MyRefTiles(:)'
                    MyBaseSpots = (ndRoundTile(:,r)==t & ndLocalTile==t2);
                    MyLocalYX = ndLocalYX(MyBaseSpots,:);
                    MyShift0 = o.RelativePos(r,:,t,t2);
                    [M, error] = PointCloudRegister(BaseYX, MyLocalYX, MyShift0, 3);
                    fprintf('Point cloud: round %d base %d tile %d->%d, error %f\n', ...
                        r, b, t2, t, error);
                    
                    MyCorrectedYX = round([MyLocalYX, ones(sum(MyBaseSpots),1)]*M);
                    ndCorrectedYX(MyBaseSpots,:,r,b) = MyCorrectedYX;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyCorrectedYX');

                end                
            else
                %BaseImSm = imdilate(BaseIm, strel('disk', 2));
                %BaseImSm = BaseIm;

                % add one because we are using coordinates starting at 0
                ndSpotColors(MySpots,b,r) = IndexArrayNan(BaseImSm, 1+ndRoundYX(MySpots,:,r)');

            end
        end
        TifObj.close();
    end
end

%% now find those that were detected in all tiles
Good = all(isfinite(ndSpotColors(:,:)),2);
GoodGlobalYX = ndGlobalYX(Good,:);
GoodSpotColors = ndSpotColors(Good,:,:);
GoodLocalTile = ndLocalTile(Good);
GoodIsolated = ndIsolated(Good);
%% plot those that were found and those that weren't
if o.Graphics
    PlotScale = 1;
    figure(1003); clf; hold on; set(gca, 'color', 'k');
    plot(ndGlobalYX(Good,2), ndGlobalYX(Good,1), 'b.', 'markersize', 1);
    plot(ndGlobalYX(~Good,2), ndGlobalYX(~Good,1), 'r.', 'markersize', 1);
    legend({'resolved', 'unresolved'});
    % now put on edges
    SquareX1 = [0, 0, o.TileSz];
    SquareY1 = [o.TileSz, 0, 0];
    SquareX2 = [o.TileSz, o.TileSz, 0];
    SquareY2 = [0, o.TileSz, o.TileSz];

    SquareColors = hsv2rgb([(1:5)'/5, [.5 0 .5 .5 .5]', [.6 1 .6 .6 .6]']);
    for r=1:2;%o.nRounds
        for t=Tiles
            plot((SquareX1 + o.RefPos(t,2)-o.RelativePos(r,2,t,t)),...
                (SquareY1 + o.RefPos(t,1)-o.RelativePos(r,1,t,t)),...
                '--', 'Color', SquareColors(r,:));
            plot((SquareX2 + o.RefPos(t,2)-o.RelativePos(r,2,t,t)),...
                (SquareY2 + o.RefPos(t,1)-o.RelativePos(r,1,t,t)),...
                ':', 'Color', SquareColors(r,:));

            text((o.RefPos(t,2)-o.RelativePos(r,2,t,t)), ...
                (o.RefPos(t,1)-o.RelativePos(r,1,t,t)), ...
                sprintf('T%d r%d', t, r), 'color', SquareColors(r,:)); 
        end
    end
    
    
    %set(gca, 'YDir', 'reverse');
end
       


%% sanity check
plsz = 7;
if o.Graphics ==2
    
    GoodRoundYX = ndRoundYX(Good,:,:);
    GoodRoundTile = ndRoundTile(Good,:);
    GoodCorrectedYX = ndCorrectedYX(Good,:,:,:);
    % which ones to display? Choose those that have different tiles on diff
    % rounds
    % PlotSpots = find(~all(bsxfun(@eq, ndRoundTile(:,1), ndRoundTile),2) & all(isfinite(ndRoundTile),2) & ndIsolated);
    % PlotSpots = 10000:10010;
    
    % roi = [y1 y2 x1 x2];
    roi = [1742 1755 213 227];
    %roi = [1342 1391 1175 1252]*o.DapiScaleFac;
    %roi = [87 107 2003 2012]*o.DapiScaleFac;
    %roi = [-inf inf -inf inf];
    PlotSpots = find(GoodGlobalYX(:,1)>roi(1) & GoodGlobalYX(:,1)<roi(2) & GoodGlobalYX(:,2)>roi(3) & GoodGlobalYX(:,2)<roi(4));
    %PlotSpots = (GoodGlobalYX(:,1)>roi(1)*4 & GoodGlobalYX(:,1)<roi(2)*4 & GoodGlobalYX(:,2)>roi(3)*4 & GoodGlobalYX(:,2)<roi(4)*4);
    %PlotSpots = PlotSpots & any(bsxfun(@ne,ndRoundTile(Good,:), ndLocalTile(Good,:)),2);
    
    for s=(PlotSpots(:))' %PlotSpots(randperm(length(PlotSpots)))'
        figure(91); clf
        for r=1:o.nRounds
            t=GoodRoundTile(s,r);

            fprintf('Spot %d, round %d, tile %d: y=%d, x=%d\n', s, r, t, GoodRoundYX(s,1,r), GoodRoundYX(s,2,r));

            Ylegends = {'Anchor', o.bpLabels{:}};
            for b=0:o.nBP
                
                      
                if b==0                    
                    y0 = GoodRoundYX(s,1,r);
                    x0 = GoodRoundYX(s,2,r);
                else
                    y0 = GoodCorrectedYX(s,1,r,b);
                    x0 = GoodCorrectedYX(s,2,r,b);
                end
                if ~isfinite(x0) | ~isfinite(y0)
                    continue;
                end
                y1 = max(1,y0 - plsz);
                y2 = min(o.TileSz,y0 + plsz);
                x1 = max(1,x0 - plsz);
                x2 = min(o.TileSz,x0 + plsz);
           
                
                BaseIm = imread(o.TileFiles{r,t}, o.AnchorChannel + b, 'PixelRegion', {[y1 y2], [x1 x2]});
                if o.SmoothSize
                    BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
                else
                    BaseImSm = BaseIm;
                end

                subplot(o.nBP+1, o.nRounds, (b)*o.nRounds + r)
                imagesc([x1 x2], [y1 y2], BaseImSm); hold on
%                imagesc([x1 x2], [y1 y2], BaseIm); hold on
                axis([x0-plsz, x0+plsz, y0-plsz, y0+plsz]);
                plot(xlim, [y0 y0], 'w'); plot([x0 x0], ylim, 'w');
                caxis([0 o.DetectionThresh*2]);
                if r==1; ylabel(Ylegends{b+1}); end
                colorbar;
                
                title(sprintf('Round %d, Base %d, Tile %d', r, b, t));
                drawnow
            end
        end
        fprintf('\n');
        figure(92); clf
        imagesc(sq(GoodSpotColors(s,:,:)));
        set(gca, 'ytick', 1:4); set(gca, 'yticklabel', o.bpLabels);
        %caxis([0 o.DetectionThresh*2]);
%         fprintf('local YX = (%f, %f) screen YX = (%f, %f) Called as %s, %s, quality %f\n', ...
%             GoodRoundYX(s,1), GoodRoundYX(s,2), GoodGlobalYX(s,1)/4, GoodGlobalYX(s,2)/4, ...
%             GoodCodes{s}, GoodGenes{s}, GoodMaxScore(s));
        pause;
    end
end



%%
save Spots GoodGlobalYX GoodSpotColors GoodIsolated