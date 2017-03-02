%% files to read
samples = ls;
samples =  cellstr(samples(3:14,:));
windowsize = 1000:1000:10000;

%% read and import
All = cell(length(samples),3);
RiskGroup = [];sy
for i = 1:length(samples)
    decodedir = strcat(samples{i});
    decodefile = strcat(decodedir, '\QT_0.4_0.005_details.csv');
    [name,pos,quality] = getinsitudata_f(decodefile);
    Risk = cell(length(windowsize),1);
    Polygon = cell(length(windowsize),1);
    
    for j = 1:length(windowsize)
        binned = binreads_f(pos, windowsize(j));
        RiskBin = cell(size(binned,2),1);
        for k = 1:length(RiskBin)
            [risk, counts] = oncotypedx_f(name(binned(:,k)));
            RiskBin{k} = risk;
        end
        Risk{j} = RiskBin;
        
        polycoord = polygonposition_f(pos, windowsize(j), 0.2);
        Polygon{j} = polycoord;
    end
    
    disp('Preparing plots..')
    imgfile = ['HE\' samples{i} '.png'];
    if exist(imgfile, 'file')
        I = imread(imgfile);
        I = imresize(I,1);
    else
        I = 0;
    end
    oncotypeplot_f(Risk, Polygon, windowsize*0.2, 0.2, I);
    mkdir('OncotypeDX\');
    saveas(gcf,['OncotypeDX\' samples{i} '.fig']);
    
%     RiskGroup = [RiskGroup; [samples{i},{risk}]];
%     All(i,:) = [{name},{pos},{quality}];
end

