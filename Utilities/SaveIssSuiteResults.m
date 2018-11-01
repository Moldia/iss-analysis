% save results from iss suite in csv format
% Xiaoyan, 2017


% load oExtra_genes_161220KI_3-1

% all reads
towrite = [o.CharCodes(o.SpotCodeNo), o.GeneNames(o.SpotCodeNo), num2cell(fliplr(o.SpotGlobalYX)), num2cell(o.SpotScore)]';

% details
fid = fopen('Basecalling_NewRegistration\SpotCoordinates.csv', 'w');
fprintf(fid, 'code,name,posX,posY,score\n');
fprintf(fid, '%s,%s,%d,%d,%d\n', towrite{:});
fclose(fid);

% count data
[~, fname] = details_to_count('Basecalling_NewRegistration\SpotCoordinates.csv', 1, 2);
movefile(fname, 'Basecalling_NewRegistration\SpotCoordinates_GeneCount.csv');

[~, fname] = details_to_count('Basecalling_NewRegistration\SpotCoordinates.csv', 1, 1);
movefile(fname, 'Basecalling_NewRegistration\SpotCoordinates_CodeCount.csv');

% details after quality threshold applied 
towrite = towrite(:,o.SpotScore>o.CombiQualThresh & o.SpotIntensity>o.CombiIntensityThresh);

fid = fopen(['Basecalling_NewRegistration\SpotCoordinates_' num2str(o.CombiQualThresh) '_' num2str(o.CombiIntensityThresh) '.csv'], 'w');
fprintf(fid, 'code,name,posX,posY,score\n');
fprintf(fid, '%s,%s,%d,%d,%d\n', towrite{:});
fclose(fid);

% gene count data after quality threshold applied
[~, fname] = details_to_count(['Basecalling_NewRegistration\SpotCoordinates_' num2str(o.CombiQualThresh) '_' num2str(o.CombiIntensityThresh) '.csv'], 1, 2);
movefile(fname, 'Basecalling_NewRegistration\GeneCount_QT.csv');
 
% code count data after quality threshold applied
[~, fname] = details_to_count(['Basecalling_NewRegistration\SpotCoordinates_' num2str(o.CombiQualThresh) '_' num2str(o.CombiIntensityThresh) '.csv'], 1, 1);
movefile(fname, 'Basecalling_NewRegistration\CodeCount_QT.csv');

