function [ccg_out, Pairs, gs, cum_dens, rel_cum_dens, ripley, local_density] = CCG_2d(v0, G0, MaxDist, nBins, gs)
% [ccg, Pairs, gs, cum_dens, rel_cum_dens, ripley] = CCG_2d(v, G, MaxDist, nBins, GSubset)
% constructs a radial 2d cross-corrlogram assuming circular symmetry
% 
% xv is a Nx2 array giving the x and y coordinates of each point
% G says which point is in which group.
% MaxDist gives the maximum distance to consider between point pairs 
% nBins: number of bins to compute the ccg
% GSubset says which groups to compute the ccg for
%
% The output array will be 3d with the first dim being radial distance and the second two 
% specifying the 2 groups in question (in the order of GSubset). This
% output is the density of a point of class two at a given radius from a
% point of class 1
%
% optional output argument pairs gives a nx2 array with the indices of the spikes
% in each train that fall in the CCG.
%
% optional output gs is a struct that returns the group labels used 
%
% optional output cum_dens gives the cumulative densitity: the density
% inside a cirlce of raidus r. rel_cum_dens gives the same thing but 
% normalized to be 1 for a homogeneous Poisson process
% optional output ripley gives Ripley's L (google it), normalized to be 1
% for a homogeneous Poisson process

if isempty(G0); G0=1; end;
if size(G0) == [1,1]
    G0 = repmat(G0, [size(v0,1) 1]);
end

% take subset 
if nargin<5; gs = unique(G0); end
[Included, G1] = ismember(G0, gs);
v = v0(Included,:);
G = G1(Included);
nG = length(gs);


% now we divide the points up into tiles of size MaxDist. We are going to
% do the computations within each tile and its neighbors to save time.

% bottom left most edge of grid. We subtract MaxDist so there is one empty tile to 
% the left of everything, so the later loops don't crash.

xMin = min(v(:,1))-MaxDist;
yMin = min(v(:,2))-MaxDist;

% compute which tile each point is in
xTile = 1+floor((v(:,1)-xMin)/MaxDist);
yTile = 1+floor((v(:,2)-yMin)/MaxDist);

% again, add empty tile to the top right
nxTiles = max(xTile)+1;
nyTiles = max(yTile)+1;

% make lists of which points in which tiles
TilePoints = cell(nxTiles, nyTiles);
for i=1:nxTiles
    for j=1:nyTiles
        TilePoints{i,j} = find(xTile==i & yTile==j);
    end
end

% now main loop. go through tiles and compute CCG of each tile with its
% neighbors. Exclude edges (which are empty by design)
local_pairs = cell(nxTiles, nyTiles);
ccg_out = zeros(nBins,nG,nG);
for i=2:nxTiles-1
    for j=2:nyTiles-1
        % find points in center tile and its neighbors
        cPoints = TilePoints{i,j};
        nPoints = unique(vertcat(TilePoints{i-1:i+1, j-1:j+1}));
        
        if isempty(nPoints) || isempty(cPoints)
            continue;
        end
        
        % now compute distance of each pair
        dx = bsxfun(@minus, v(cPoints,1), v(nPoints,1)');
        dy = bsxfun(@minus, v(cPoints,2), v(nPoints,2)');
        r = sqrt(dx.^2 + dy.^2);
        r(bsxfun(@eq, cPoints(:), nPoints(:)')) = Inf; % don't count selfies
        
        % find neighbors that are close enough
        [lc, ln] = find(r<MaxDist);
        local_pairs{i,j} = [cPoints(lc,1), nPoints(ln,1)]; %,1 to make it a column
        
        % make histogram of distances for each close pair
        local_r = r(sub2ind(size(r), lc, ln));
        local_gps = G(local_pairs{i,j});
        if size(local_gps,2)==1; local_gps=local_gps'; end; % GOD I HATE MATLAB
        % ccg is indexed by r, then two group labels
        ccg_subs = [1+floor(local_r/MaxDist*nBins)];
        ccg_out = ccg_out + accumarray([ccg_subs(:), local_gps],1, [nBins nG nG]);
    end
end

% now we want to divide the ccg by the area of each annulus
area = ((1:nBins).^2 - (0:nBins-1).^2)'*pi*(MaxDist/nBins).^2;
ccg_out_norm1 = bsxfun(@rdivide, ccg_out, area);

% now we want to divide by the number of points of each class
n_in_class = accumarray(G,1);
local_density = bsxfun(@rdivide, ccg_out_norm1, n_in_class(:)');

% now put together all pairs - and convert back to original numbering
% (including groups that were not computed)
all_pairs = vertcat(local_pairs{:});
reindex_array = find(Included);
Pairs = reindex_array(all_pairs);
        
if nargout>=4
    Total_Area = prod(max(v0) - min(v0));
    CumDist = cumsum(ccg_out, 1);
    Norm1 = bsxfun(@rdivide, CumDist, reshape(n_in_class, [1, nG, 1]));
    Norm2 = bsxfun(@rdivide, Norm1, reshape(n_in_class/Total_Area, [1, 1, nG]));
    r = (1:nBins)'*MaxDist/nBins;
    rel_cum_dens = bsxfun(@rdivide, Norm2, pi*r.^2);
    cum_dens = bsxfun(@rdivide, Norm1, pi*r.^2);
    ripley = bsxfun(@rdivide, sqrt(Norm2/pi*Total_Area), r);
end
    