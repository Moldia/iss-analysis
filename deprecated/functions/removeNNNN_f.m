function [name_uni, idx_name, pos] = removeNNNN_f(name, pos)

[name_uni, ~, idx_name] = unique(name);
idxNNNN = strcmp(name_uni,'NNNN');

if nnz(idxNNNN)
    idxNNNN = find(idxNNNN);
    pos = pos(idx_name~=idxNNNN,:);
    name = name(idx_name~=idxNNNN);
    [name_uni, ~, idx_name] = unique(name);
end

end

