usedbarcodes = importdata('C:\Users\Xiaoyan\SkyDrive\usedbarcodes.csv');
usedbarcodes = cellfun(@(v) strsplit(v, ','), usedbarcodes, 'uni', 0);
usedbarcodes = reshape([usedbarcodes{:}], 2, [])';
usedbarcodes = usedbarcodes(:,2);
usedbarcodes = unique(usedbarcodes);

lenbarcode = cellfun(@length, usedbarcodes);

DO = cellfun(@(v) str2num(v(1)), usedbarcodes(lenbarcode==5), 'uni', 0);
DOorder = {'T' 'A' 'C' 'G' 'N' 'N'};
DO = DOorder(cell2mat(DO));
barcodes = cellfun(@(v) v(2:end), usedbarcodes(lenbarcode==5), 'uni', 0);
barcodes = strcat(barcodes, DO');

[~, letterset] = barcode_generator(5,4,barcodes,2);

newbarcodes = letterset{2};
newbarcodes = cellfun(@(v) [num2str(find(strcmp(DOorder, v(5)))), v(1:4)],...
    newbarcodes, 'uni', 0);
