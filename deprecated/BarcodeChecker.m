barcodes = importdata('C:\Users\qxyyx\OneDrive\worktemp\Probe_all_temp.xlsx');
barcodes = barcodes.textdata.Sheet1;


%% with DO
% DO = barcodes(:,1);
% DO = cellfun(@(v) str2num(v(3)),DO);

barcodes = barcodes(:,6);
% remove header
barcodes = barcodes(2:end);
DO = cellfun(@(v) str2num(v(1)),barcodes,'uni',0);

normalplp = cellfun(@(v) ~isempty(v), DO);
DO = cell2mat(DO);
barcodes = cellfun(@(v) v(2:5),barcodes(normalplp),'uni',0);
barcodes = letter2num(barcodes,4);

barcodes = [DO, barcodes(:,2:end)];

%% no DO
% barcodes = letter2num(barcodes,4);
% barcodes = barcodes(:,2:end);

%%
Dist = zeros(length(barcodes),'uint8');
for i = 1:5
    Dist = Dist + uint8((repmat(barcodes(:,i),1,length(barcodes)) ~= repmat(barcodes(:,i)',length(barcodes),1)));
end
Dist(Dist>2) = 2;

figure,imshow(Dist,[]);
colormap hot
axis on