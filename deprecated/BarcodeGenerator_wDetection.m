

%% generate list of N base barcodes
num_hybs = 4;
num_bases = 4;
[Y{num_hybs:-1:1}] = ndgrid(1:num_bases);
list = reshape(cat(num_hybs+1,Y{:}),[],num_hybs);

letters = num2letters_2(list);
clear Y
Letters_set = {};


%% homopolymers
homomer = {'AAAA' 'CCCC' 'GGGG' 'TTTT'};
homomer = cellfun(@(v) strcmp(v,cellstr(letters)),homomer,'uni',0);
homomer = cellfun(@(v) find(v),homomer);
homomer = accumarray(homomer',1,[num_bases^num_hybs 1]);


%% detection sites
num_DO = 6;
list_DO = repmat(1:num_DO,num_bases^num_hybs,1);
list_DO = list_DO(:);
list = [list_DO,repmat(list,num_DO,1)];

letters = [num2str(list_DO),repmat(letters,num_DO,1)];

homomer = repmat(homomer,num_DO,1);


%% already used in current scheme
occupied = importdata('barcodes_used.xlsx');
occupied = occupied.Sheet1(:,3);

occupied = cellfun(@(v) strcmp(v,cellstr(letters)),occupied,'uni',0);
occupied = cellfun(@(v) find(v),occupied,'uni',0);
occupied = cell2mat(occupied); % remove special barcodes, e.g. sandwich barcode
occupied = logical(accumarray(occupied,1,[num_bases^num_hybs*num_DO 1]));

%% calculate all pair-wise distances
base1 = repmat(list(:,1),1,num_bases^num_hybs*num_DO) ~= repmat(list(:,1)',num_bases^num_hybs*num_DO,1);
base2 = repmat(list(:,2),1,num_bases^num_hybs*num_DO) ~= repmat(list(:,2)',num_bases^num_hybs*num_DO,1);
base3 = repmat(list(:,3),1,num_bases^num_hybs*num_DO) ~= repmat(list(:,3)',num_bases^num_hybs*num_DO,1);
base4 = repmat(list(:,4),1,num_bases^num_hybs*num_DO) ~= repmat(list(:,4)',num_bases^num_hybs*num_DO,1);
base5 = repmat(list(:,5),1,num_bases^num_hybs*num_DO) ~= repmat(list(:,5)',num_bases^num_hybs*num_DO,1);

Dist = base1 + base2 + base3 + base4 + base5;

% figure;
% imshow(Dist,[])

%% at least two-base difference between each pair
% and at least two-base difference from current barcode list
Letters_set = [];
occupied_neighbor = Dist(occupied,:);
occupied_neighbor = sum(occupied_neighbor==1,1);
occupied_neighbor = accumarray((find(occupied_neighbor))',1,[num_bases^num_hybs*num_DO 1]);

Final = 1:2;    % in order to escape loop earlier if only fragmented barcode set left
while nnz(~occupied) && length(Final)>1
    Discard = false(num_bases^num_hybs*num_DO,1);
    Discard = Discard | occupied | homomer | occupied_neighbor;
    for i = 1:num_bases^num_hybs*num_DO
        if ~Discard(i)
            temp = Dist(i,:)==1;
            Discard = Discard | temp';
        end
    end
    
    
    occupied = occupied | ~Discard;
    Final = find(~Discard);
    Letters_set = [Letters_set,{cellstr(letters(Final,:))}];
end
