function show_only_legend(ah, legendnames, symbols)
% show_only_legend(ah, legendnames, symbols)
% create an invisible dummy figure in subplot to show only legend
% Xiaoyan, 2017


randomPoints = rand(length(legendnames),2);

axes(ah);
hold on;
for i = 1:numel(legendnames)
    ph = plot(randomPoints(i,1), randomPoints(i,2), symbols{i});
    set(ph, 'visible', 'off');
end
h = legend(legendnames, 'color', [.6 .6 .6], 'fontsize', 5);
set(h, 'location', 'west');
set(ah, 'visible', 'off');
end
