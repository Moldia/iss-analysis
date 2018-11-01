% extract read coordinates from a .fig file
% Xiaoyan, 2017

scale = .25;    % current scale
ch = get(gca, 'child');
total = 0;
for i = 1:length(ch)
    if strcmp(ch(i).Type, 'line')
        total = [total; total(end)+length(ch(i).XData)];
    end
end

pos = zeros(total(end),2);
nameidx = zeros(total(end),1);
names = [];
for i = 1:length(total)-1
    pos(total(i)+1:total(i+1),:) = [ch(i).XData', ch(i).YData'];
    nameidx(total(i)+1:total(i+1)) = i;
    names = [names; {ch(i).DisplayName}];
end
pos = pos/scale;

fid = fopen('SpotCoordinates.csv', 'w');
fprintf(fid, 'Name,posX,posY\n');
towrite = [names(nameidx), num2cell(pos)]';
fprintf(fid, '%s,%d,%d\n', towrite{:});
fclose(fid);

I = ch(end).CData;
imwrite(I, ['DAPI_', num2str(scale*100), '%.tif']);