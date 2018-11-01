%% get the list of barcodes
%  with SW probes
%  Xiaoyan, 2017

%% load database file and used probe list file
allProbes = importdata('C:\Users\Xiaoyan\OneDrive\worky\allProbes.xlsm');
allProbes = allProbes.textdata.cDNAprobes;
allProbes = allProbes(2:end, :);

qProbes = importdata('probesused.csv');
seqorder = '12345';
baseChannel = [3,4,5,6];      % channel order of A, C, G, T
% DObase = {'T' 'A' 'C' 'G' 'N' 'N'};    % channel order of detection oligos (DO1-DO6)
nSWChannels = 6;

[Codebook, plProbes, barcodeLetter, swProbes] = ...
    get_codebook(qProbes, allProbes(:,1), allProbes(:,5), allProbes(:,8),...
    seqorder, baseChannel);

%% write files
fid = fopen('IDlist.csv', 'w');

barcodeDigit = [allProbes(plProbes,5), barcodeLetter];

[~, idxUni] = unique(barcodeDigit(:,2), 'stable');
barcodeDigit = (barcodeDigit(idxUni,[2,1]))';
fprintf(fid, '%s,%s\n', barcodeDigit{:});
fclose(fid);

barcodeAll = [barcodeLetter; allProbes(swProbes,6)];
probesAll = [plProbes; swProbes];
probesAll = allProbes(probesAll, [1,4,5,6]);

% unique barcodes
[~, idxUni] = unique(barcodeAll, 'stable');
fid = fopen('codebook_unique.csv', 'w');
fprintf(fid, 'probeID,name,gene,barcode,barcode_bases');
for i = 1:length(seqorder)
    for j = 1:6
        fprintf(fid, strcat(',base', num2str(i), '_channel', num2str(j)));
    end
end
for j = 1:nSWChannels
    fprintf(fid, strcat(',SW_channel', num2str(j)));
end
fprintf(fid, '\n');
fmt = ['%s,%s,%s,%s,%s,', lineformat('%d', 6*length(seqorder)+nSWChannels)];
barcodeDigit = ([probesAll(idxUni,:), barcodeAll(idxUni), num2cell(Codebook(idxUni,:))])';
fprintf(fid, fmt, barcodeDigit{:});
fclose(fid);

% fid = fopen('codebook_all.csv', 'w');
% fprintf(fid, 'probeID,name,gene,barcode,barcode_bases');
% for i = 1:5
%     for j = 1:6
%         fprintf(fid, strcat(',base', num2str(i), '_channel', num2str(j)));
%     end
% end
% for j = 1:nChannelsSW
%     fprintf(fid, strcat(',SW_channel', num2str(j)));
% end
% fprintf(fid, '\n');
% fmt = ['%s,%s,%s,%s,%s,', lineformat('%d', 6*length(seqorder)+nChannelsSW)];
% barcodeDigit = ([probesAll, barcodeAll, num2cell(Codebook)])';
% fprintf(fid, fmt, barcodeDigit{:});
% fclose(fid);

