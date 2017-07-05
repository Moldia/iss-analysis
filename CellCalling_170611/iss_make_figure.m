function iss_make_figure(o, SpotsYX, Genes, Sizes)
% InSituPlot(o, SpotsYX, Genes, Sizes)
%
% plot the results of in situ sequencing spot detection. 
% SpotYX gives coordinates; Gene is a list of strings; 
% 
% BackgroundImage is loaded using options o
%
% sizes can be a vector or a scalar - only used for scatter, which isn't
% called anyway.
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
if nargin<4; Sizes = 6; end;

Scatter=0; % faster

clf; set(gcf, 'color', 'k');
set(gca, 'color', 'k');
hold on
% 
% if nargin>=3
%     if isstr(BackgroundImage)
%         BackgroundImage = imread(BackgroundImage);
%     end
%     imagesc(BackgroundImage);
% end
if exist(fullfile(o.OutputDirectory, 'background_boundaries.tif'), 'file')==2
    BackgroundImage = imread(fullfile(o.OutputDirectory, 'background_boundaries.tif'));
else
    BackgroundImage = imread(fullfile(o.OutputDirectory, 'background_image.tif'));
end

imagesc(BackgroundImage);
colormap(bone);

uGenes = unique(Genes);

if Scatter
    for i=1:length(uGenes)
        g = uGenes(i);
        MySpots = strcmp(Genes, g);
        h(i) = scatter(SpotsYX(MySpots,2)/o.DapiScaleFac, SpotsYX(MySpots,1)/o.DapiScaleFac, Sizes(MySpots));
    end
else
    for i=1:length(uGenes)
        g = uGenes(i);
        MySpots = strcmp(Genes, g);
        h(i) = plot(SpotsYX(MySpots,2)/o.DapiScaleFac, SpotsYX(MySpots,1)/o.DapiScaleFac, '.');
    end 
end

hl = legend(h,uGenes);

% gscatter(SpotsYX(:,2), SpotsYX(:,1), Genes, [], [], Sizes)
iss_change_gene_symbols(0,[],2);

end

