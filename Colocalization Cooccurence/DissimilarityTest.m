% test the null hypothesis if two groups of data points are drawn from the same continuous distribution 
% two-sample Kolmogorov-Smirnov test
% Xiaoyan, 2017

[name, pos] = getinsitudata('Decoding\QT_0.5_0_details.csv');
[uniName, ~, idxName] = unique(name);
histName = hist(idxName, 1:max(idxName));

% remove NNNN and any reads<30
[name, pos] = removereads(name, [uniName(histName<30); {'NNNN'}], pos);
[uniName, ~, idxName] = unique(name);

% bin
pos_bin = ceil(pos/200);
pos_bin = max(pos_bin(:,1))*(pos_bin(:,1)-1) + pos_bin(:,2);


ksstat = zeros(length(uniName), 2);
for i = 1:length(uniName)
    for j = i:length(uniName)
        [~, p, ks2stat] = kstest2(pos_bin(idxName==i), pos_bin(idxName==j));
        ksstat(i,j,1) = ks2stat;
        ksstat(i,j,2) = p;
        ksstat(j,i,1) = ks2stat;
        ksstat(j,i,2) = p;
    end
end
figure; imagesc(ksstat(:,:,2));
set(gca, 'xtick', 1:length(uniName), 'xticklabel', uniName,...
    'ytick', 1:length(uniName), 'yticklabel', uniName,...
    'xticklabelrotation', 90);


% sorting
L = linkage(ksstat(:,:,2));
order = dendroperm(L);

figure; imagesc(ksstat(order,order,2));
set(gca, 'xtick', 1:length(uniName), 'xticklabel', uniName(order),...
    'ytick', 1:length(uniName), 'yticklabel', uniName(order),...
    'xticklabelrotation', 90);

c = colorbar;
c.Label.String = 'Kolmogorov-Smirnov test p-value';
c.Label.FontSize = 12;