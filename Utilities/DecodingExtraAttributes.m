%% add extra column based on other infomration
%  Xiaoyan, 2016-5-25

[name,pos,tileid]= getinsitudata_f('QT_0.4_0.001_details.csv',2,1,3);
% [code,genral_quality]= getinsitudata_f('QT_0.4_0.001_details.csv',1,4);
IHC_bw = imread('IHC_epithelial_final_20%.png');

pos_scaled = correctcoordinates_f(pos,0.2);
pos_pix = round(pos_scaled);
pos_linear = (pos_pix(:,1)-1)*size(IHC_bw,1)+pos_pix(:,2);
tumor = IHC_bw(pos_linear)==255;

towrite = [name,num2cell([pos,double(tumor)])];

% towrite = [code(tumor),name(tumor),num2cell([pos(tumor,:),tileid(tumor),genral_quality(tumor,:)])];



fid = fopen('K_QT_wTumor.csv','w');
fprintf(fid,'name,pos_x,pos_y,tumor\n');
towrite = towrite';
fprintf(fid,'%s,%d,%d,%d\n',towrite{:});
fclose(fid);