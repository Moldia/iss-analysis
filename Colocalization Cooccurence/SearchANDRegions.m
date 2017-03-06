% search regions where all specified genes co-occur
% Xiaoyan, 2017

[name, pos] = getinsitudata('C:\Users\qxyyx\Desktop\temp\beforeQT_details.csv');

cooccur = search_reads_cooccur(name, pos, 100);
cooccur = search_reads_cooccur(name, pos, 100, {'Drd3', 'Prrx1', 'Ctdnep1'});


