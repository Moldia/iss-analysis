% compare plot and Allen ISH images side by side
% Xiaoyan, 2017


% number of probes per gene
libProbes = importdata('C:\Users\Xiaoyan\OneDrive\worky\allProbes.xlsm');
allProbes = libProbes.textdata.Sheet1;
allProbes = allProbes(2:end,[1,5]);
[uTargets, nProbes] = find_oligo_genes('..\probesused.csv', allProbes);


[name, pos] = getinsitudata('..\iss\SpotCoordinates.csv', 1);
[uNames, ~, idxName] = unique(name);

I = [];
for i = 1:numel(uNames)
    try
        files = cellstr(ls(fullfile('E:\PROOOJECTS\12_Neuron_mapping\Genes\AllenISH\LowResCA1', ...
            uNames{i})));
        ISH = imread(fullfile('E:\PROOOJECTS\12_Neuron_mapping\Genes\AllenISH\LowResCA1', ...
            uNames{i}, files{3}));
        
        idx = find(strcmp(uNames{i}, uTargets));
        
        clf;
        subplot(121)
        plot(pos(idxName==i,1), pos(idxName==i,2), '.', 'markersize', 1);
        title(['nProbes = ' num2str(nProbes(idx))])
        axis image
        axis off
        
        subplot(122)
        imshow(ISH)
        title(uNames{i})
        
        drawnow
        I = [I, {frame2im(getframe(1))}];
        saveas(gcf, [uNames{i} '.jpg']);
    end
end

% gif
for i = 1:numel(I)
    [A, map] = rgb2ind(I{i}, 256);    
    if i == 1
        imwrite(A, map, 'AllenISH_170315_161220KI_4-1.gif', 'gif', 'Loopcount', inf, 'DelayTime', .5);
    else
        imwrite(A, map, 'AllenISH_170315_161220KI_4-1.gif', 'gif', 'WriteMode', 'append', 'DelayTime', .5);
    end
end

