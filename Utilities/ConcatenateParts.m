%% sequencing
for b = 1:4
    b
    for c = 1:5
        c
        I1 = imread(['base' num2str(b) '_upper\' 'base' num2str(b) '_upper_c' num2str(c) '_ORG.tif']);
        I2 = imread(['base' num2str(b) '_mid\' 'base' num2str(b) '_mid_c' num2str(c) '_ORG.tif']);
        I3 = imread(['base' num2str(b) '_lower\' 'base' num2str(b) '_lower_c' num2str(c) '_ORG.tif']);
        I = [I1;I2;I3];
        imwrite(I,['base' num2str(b) '\base' num2str(b) '_c' num2str(c) 'fullsize.tif']);
    end
end

% %% single image
% for c = 1:3
%     I1 = imread(['HE_upper\HE_upper_c' num2str(c) '_ORG.tif']);
%     I2 = imread(['HE_mid\HE_mid_c' num2str(c) '_ORG.tif']);
%     I3 = imread(['HE_lower\HE_lower_c' num2str(c) '_ORG.tif']);
%     I = [I1;I2;I3];
%     imwrite(I,['HE\HE_c' num2str(c) '_fullsize.tif']);
% end