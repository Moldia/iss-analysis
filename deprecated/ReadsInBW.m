%% read xlsx file
samples = importdata('..\PDGF progress_170122.xlsx');
scale = samples.data(:,21);
qt = samples.data(:,20);
group = samples.textdata(2:end,2);
group = cellfun(@(v) [upper(v(1)), v(2:end)], group, 'uni', 0);
samples = strcat(samples.textdata(2:end,3), '_', samples.textdata(2:end,2));
bwtype = {'cancer' 'stroma'};
sym = {'g.' 'y.'};

%% loop through all images
imgfls = cellstr(ls('..\IHC\cancer_stroma_regions\'));
for i = 1:length(samples)
    if ~isnan(scale(i))
        backfls = cellstr(ls('..\IHC\aligned\'));
        backimg = cellfun(@(v) (~isempty(strfind(lower(v),lower(samples{i})))&~(strcmp(v(1),'.'))), backfls);
        backimg = backfls(backimg);
        background = imread(['..\IHC\aligned\' backimg{1}]);            

        % load data
        posfolder = ['..\Sequencing\' group{i} '\' samples{i}];
        posfls = cellstr(ls(posfolder));
        posfile = cellfun(@(v) (~isempty(strfind(v,'CP_'))|~isempty(strfind(v,'Cellprofiler'))), posfls);
        posfolder = strcat(posfolder, '\', posfls(posfile), '\Decoding');
        
%         posfls = cellstr(ls(posfolder{1}));
%         posfile = cellfun(@(v) (isempty(strfind(v,'beforeQT'))&~isempty(strfind(v,'_details'))), posfls);
%         posfile = posfls(posfile);
%         data = importdata(strcat(posfolder{1}, '\', posfile{1}));
        try
            data = importdata([posfolder{1} '\QT_' num2str(qt(i)) '_0.0003_details.csv']);
            [name, pos] = getinsitudata_f([posfolder{1} '\QT_' num2str(qt(i)) '_0.0003_details.csv']);
        catch
            posfolder = {['..\Sequencing\' group{i} '\' samples{i} '\Decoding']};
            data = importdata([posfolder{1} '\QT_' num2str(qt(i)) '_0.0003_details.csv']);
            [name, pos] = getinsitudata_f([posfolder{1} '\QT_' num2str(qt(i)) '_0.0003_details.csv']);
        end

        header = data.textdata(1,:);
        textdata = data.textdata(2:end,:);
                
        clf;
        ax = subplot(1,3,1); imshow(background); hold on; title(strrep(samples{i},'_','\_'));
        Ax = ax;
        for j = 1:2
            
            % find binary image
            imgfile = cellfun(@(v) (~isempty(strfind(v,bwtype{j}))&~isempty(strfind(lower(v),lower(samples{i})))), imgfls);
            imgfile = imgfls(imgfile);
            
            % load binary image
            IHC_bw = imread(['..\IHC\cancer_stroma_regions\' imgfile{1}]);
            IHC_bw = logical(IHC_bw(:,:,1));

                        
            % filter reads
            pos_temp = correctcoordinates_f(pos,scale(i));
            pos_pix = round(pos_temp);
            pos_linear = (pos_pix(:,1)-1)*size(IHC_bw,1)+pos_pix(:,2);
            specific = IHC_bw(pos_linear)==1;
           
            towrite = [textdata(specific,1:2),num2cell(data.data(specific,:))]';
            fid = fopen(strcat('ReadsInBW\', samples{i}, '_', bwtype{j}, '_details.csv'), 'w');
            fmt = repmat('%s,',1,length(header));
            fmt = [fmt(1:end-1), '\n'];
            fprintf(fid,fmt,header{:});
            fmt = ['%s,%s',repmat(',%d',1,length(header)-2),'\n'];
            fprintf(fid,fmt,towrite{:});
            fclose(fid);
            
            % visualization
            subplot(1,3,1); plot(pos_temp(specific,1),pos_temp(specific,2),sym{j});           
            ax = subplot(1,3,j+1); imshow(IHC_bw); hold on; plot(pos_temp(specific,1),pos_temp(specific,2),sym{j},'linewidth',1.5);
            title(bwtype{j});
            drawnow;
            Ax = [Ax, ax];
          
        end
        
        linkaxes(Ax, 'xy');
        saveas(gcf,['ReadsInBW\' samples{i} '.png']);
    end
end
