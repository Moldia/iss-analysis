function [L1,NN_num] = listself(position,idx,uniname_re,radius)
% used for Cluster_Area
% connects neighboring transcripts

L1 ={};

pos_x = position(uniname_re==idx,1);
pos_y = position(uniname_re==idx,2);

[I,~] = rangesearch([pos_x,pos_y],[pos_x,pos_y],2*radius);
NN_num = zeros(length(I),1);
for i = 1:length(I)
    NN_num(i) = size(I{i},2);
end

I2 = I;
I2(NN_num==1) = {[]};
NN_num(NN_num==1) = 0;

i = 1;
while i <= length(NN_num)
    if (NN_num(i))~=0
        L = [];
        L1 = [L1;{childlist(i,L)}];
    else
        L1 = [L1;{[]}];
    end
    i = i+1;
end


% RECURSION
    function L = childlist(M,L)
        for j = 1:length(M)
                N = I2{M(j)};
                N = unique(N,'stable');
                N = N(~ismember(N,L));
                I2(M(j)) = {[]};
                L = [L,N];
                if isempty(N)
                else
                    L = childlist(N,L);
                end
        end
    end



end