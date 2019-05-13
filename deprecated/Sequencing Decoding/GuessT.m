function [data,columnn,rown]=GuessT(expect_matrix,list,name_tag_list,exp_tags)

% Make best guess when T channel is removed from analysis.
% Only one mismatch is allowed and that mismatch has to be T.
% Currently only supports one T in barcode.
%
% GuessT v2.1
% May 3, 2014, Xiaoyan


% extract unexpected reads, and make a digit matrix of unique reads
%----------------------------
NNNN = list(strcmp(name_tag_list,'NNNN'));
[numN,uni_N] = hist(NNNN,unique(NNNN));

uni_N_str = num2str(uni_N);
N_matrix = [];
for j = 1:length(uni_N_str)
    N_digits = [];
    for i = 1:length(expect_matrix(1,:))
        N_digits = [N_digits str2num(uni_N_str(j,i))];
    end
    N_matrix = [N_matrix; N_digits];
end

% "sequence alignment"
%----------------------------
for kk = 1: length(N_matrix(:,1))
    N_read = N_matrix(kk,:);
    for j = 1:length(expect_matrix(:,1))
        exp_read = expect_matrix(j,:);
        x = 0; y=0;
        for i = 1:length(expect_matrix(1,:))
            if N_read(i) == exp_read(i)
                x = x+1;
            elseif exp_read(i) == 4
                y = y+1;
            end
        end
        Ascore(kk,j) = x;
        Tscore(kk,j) = y;
    end
end
% find one mismatch and that is T
%----------------------------
f = (Ascore(:,:)==length(expect_matrix(1,:))-1 & Tscore(:,:)==1);
show=cell(length(uni_N),length(expect_matrix(:,1)));
show(f)={'match'};
show(~f)={''};
fi=figure;
for i = 1:length(numN)
    rown{i}= num2letters(uni_N(i));
end
columnwide=[];
for i = 1:length(expect_matrix(:,1))
    columnwide = [columnwide, num2cell(80)];
end
columnn = ['count'; exp_tags];
data = [num2cell(numN(:)) show];
h = uitable(fi,'data',data,'ColumnName',columnn,'RowName',rown);
set(h,'units','normalized','ColumnWidth',columnwide,'Position',[0 0 1 1]);
set(fi,'name','guess the T');

end