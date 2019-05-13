function [reads,names,counts] = getcountdata_f(countfile)

% extract data from gene_n_count or code_n_count files
% Xiaoyan 2015-8-12

countdata = importdata(countfile,',');
countdata(1,:) = [];
countdata = cellfun(@(v) strsplit(v,','),countdata,'UniformOutput',false);
reads = cellfun(@(v) v{:,1},countdata,'UniformOutput',false);
counts = cellfun(@(v) v{:,2},countdata,'UniformOutput',false);
counts = cell2mat(cellfun(@str2double,counts,'UniformOutput',false));
names = cellfun(@(v) v{:,3},countdata,'UniformOutput',false);


end
