% search regions where all specified genes co-occur
% Xiaoyan, 2017

[name, pos] = getinsitudata('Decoding\beforeQT_details.csv');

search_reads_cooccur(name, pos, 100, {'Drd3', 'Prrx1', 'Ctdnep1'});

