function taglist = formattaglist_f(taglist)
% format taglist
% Xiaoyan 2014-11-26


taglist_temp = cellfun(@strsplit,taglist(:,1),'uniformoutput',false);
sym_temp = taglist(:,2);
taglist = cell(length(taglist_temp),3);
for i = 1:length(taglist_temp) 
    taglist(i,:) = [taglist_temp{i},sym_temp(i)];
end