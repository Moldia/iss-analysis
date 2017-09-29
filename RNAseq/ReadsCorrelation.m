%% different lanes

 




FPKM = cell(4,4);





for s = 25:28
    for l = 1:4
        data = importdata(fullfile('G:\RNAseq\170825_NB501365_0108_AHJ2HNBGX3\UngappedMapping\fpkm',...
            ['S' num2str(s) '_L' paddigits(l, 3)],...
            'fpkm_sipmle.csv'));
        
        name = data.textdata(:,2);
        fpkm = data.data;
        
        FPKM{s-24,l} = fpkm;
        
        
    end
end


        