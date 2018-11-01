%% files to read
samples = ls();
files = cellstr(files(3:end,:));
windowsize = 1000:1000:10000;

%% read and import
Genes = {'ACTB','GAPDH','RPLP0','GUS1','TFRC','Ki-67','STK15','BIRC5','CCNB1','MYBL2','MMP11','CTSL2','GRB7','HER2','ER','PR','BCL2','SCUBE2','GSTM1','CD68','BAG1'};
All = cell(length(files),3);
CountsAll = cell(length(files),length(windowsize));
Polygon = cell(length(files),length(windowsize));

for i = 1:length(files)
    decodedir = strcat(files{i});
    decodefile = strcat(decodedir, '\QT_0.35_0.004_details.csv');
    [name, pos] = getinsitudata(files{i});

    % take only ODX genes
    pos = pos(ismember(name,Genes),:);
    name = name(ismember(name,Genes));
    
    % reorder based on the origianl gene list
    [uNames, ~, iName] = unique(name);
    idx_list = cellfun(@(v) strcmp(v,Genes), uNames, 'uni', 0);
    idx_list = cellfun(@find, idx_list, 'uni', 0);
    idx_list = cellfun(@(v) [v, zeros(1,1-length(v))], idx_list);
   
    % bin reads
    for j = 1:length(windowsize)
        binned = binreads(pos, windowsize(j));
        counttemp = double(binned).*repmat(iName,1,size(binned,2));
        
        count = zeros(21,size(binned,2));
        counttemp = hist(counttemp,0:length(uNames));
        counttemp = counttemp(2:end,:);
        count(idx_list(idx_list~=0),:) = counttemp;
        
        CountsAll{i,j} = count;

        polycoord = polygonposition(pos, windowsize(j), 0.2);
        Polygon{i,j} = polycoord;
    end
        
end

%% write file for R
Geneswrite = {'ACTB','GAPDH','RPLPO','GUS','TFRC','Ki67','STK15','BIRC5','CCNB1','MYBL2','MMP11','CTSL2','GRB7','HER2','ER','PGR','BCL2','SCUBE2','GSTM1','CD68','BAG1'};
for i = 1:length(windowsize)
    fid = fopen(['windowsize_' num2str(windowsize(i)) '.csv'], 'w');
 
    countsum = [];
    header = [];
    for j = 1:length(files)
        header = [header,...
            strcat(repmat(strrep(files(j), '\','/'), 1, size(CountsAll{j,i},2)), '_',...
            cellfun(@(v) num2str(v), num2cell(1:size(CountsAll{j,i},2)),'uni',0))];
        countsum = [countsum, CountsAll{j,i}];
    end
    
    % threshold: at least 20 reads in one grid
    header = header(sum(countsum,1)>=20);
    countsum = countsum(:,sum(countsum,1)>=20);

    header = [{'GeneName'}, header];
    header = strcat(header, ',');
    fmt = repmat('%d,',1,length(header)-1);
    header = [header{:}];
    header = [header(1:end-1), '\n'];
    fmt = ['%s,' fmt(1:end-1), '\n'];
    
    fprintf(fid, header);
    towrite = [Geneswrite', num2cell(countsum)]';
    fprintf(fid, fmt, towrite{:});
    fclose(fid);
end

%% call R for OncotypeDX scoring 
dir = cd();
dir = strrep(dir, '\', '/');
outdir = [dir, '/ODX_R'];
mkdir(outdir)
for i = length(windowsize):-1:1
    infile = [dir '/' 'windowsize_' num2str(windowsize(i)) '.csv'];
    outfile = [outdir, '/', 'windowsize_' num2str(windowsize(i)), '_ODX.xlsx'];
    rcommand = ['Rscript --vanilla E:/GitLocal/iss-analysis/Subtyping/ODX_Xiaoyan.R ',...
        infile ' ' outfile];
    system(rcommand);
end
    
%% read ODX classification from R output files
RiskGroup = cell(length(files),length(windowsize));

for i = 2:length(windowsize)   
    outfile = [outdir, '/', 'windowsize_' num2str(windowsize(i)), '_ODX.xlsx']; 
    [~, rsdata] = xlsread(outfile, 1);
    gridid = rsdata(2:end,1);
    gridid = cellfun(@(v) v(2:end), gridid, 'uni', 0);
    gridid = cellfun(@(v) strsplit(v,'_'), gridid, 'uni', 0);
    
    for j = 1:length(files)
        temp = cellfun(@(v) strcmp([v{1} '_' v{2}],files{j}), gridid);
        tempgrid = [gridid{temp}];
        tempgrid = reshape(tempgrid,3,[]);
        tempgrid = cellfun(@str2num, tempgrid(3,:));
        Temp = repmat({'NA'},max(tempgrid),1);
        Temp(tempgrid) = rsdata(temp,end);
        RiskGroup{j,i} =  Temp;
    end
end

%% figures
for i = 1:length(files)
    disp('Preparing plots..')
    imgfile = ['HE\' files{i} '.png'];
    if exist(imgfile, 'file')
        I = imread(imgfile);
        I = imresize(I,1);
    else
        I = 0;
    end
    oncotypeplot_f(RiskGroup(i,2:end), Polygon(i,2:end), windowsize(2:end)*0.2, 0.2, I);
    mkdir('OncotypeDX\');
    saveas(gcf,['OncotypeDX\' files{i} '.fig']);
    
%     RiskGroup = [RiskGroup; [samples{i},{risk}]];
%     All(i,:) = [{name},{pos},{quality}];
end

