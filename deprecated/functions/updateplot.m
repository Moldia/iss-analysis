function updateplot(f, list, varargin)

ax = get(f,'currentaxes');
h = ax.Children;

H = [];
for i = 1:length(h)
    if strcmp(class(h(i)), 'matlab.graphics.chart.primitive.Line')
        H{i} = h(i);
    end
end

if isempty(varargin)
    varargin = 'gene';
end

col_l = {'r' 'g' 'y' 'm' 'c' 'w' 'b'};
col_rgb = {'1  0  0' '0  1  0' '1  1  0' '1  0  1' '0  1  1' '1  1  1' '0  0  1'};

switch varargin
    case 'gene'
        if strcmp(list, 'all')
            for i = 1:length(H)
                H{i}.Visible = 'on';
            end
        else  
            for i = 1:length(H)
                H{i}.Visible = 'off';
            end
            G = cellfun(@(v) strcmp(v.DisplayName, list), H, 'uni', 0);
            G = cell2mat(G');
            for i = 1:length(list)
                t = find(G(:,i)==1);
                if ~isempty(t)
                    H{t}.Visible = 'on';
                    fprintf('%10s\t%s %s\n',...
                        H{t}.DisplayName, col_l{strcmp(num2str(H{t}.Color),col_rgb)}, H{t}.Marker);
                else
                    warning([list{i} ' not found']);
                end
            end
        end

        hLegend = findall(f,'tag','legend');
        uic = get(hLegend,'UIContextMenu');
        uimenu_refresh = findall(uic,'Label','Refresh');
        callback = get(uimenu_refresh,'Callback');
        hgfeval(callback,[],[]);
        
    case 'symbol'
        if strcmp(list, 'all')
            for i = 1:length(H)
                H{i}.Visible = 'on';
            end
        else
            sym = cellfun(@(v) v(2), list, 'uni', 0);

            for i = 1:length(H)
                H{i}.Visible = 'off';
            end
            G = cellfun(@(v) strcmp(v.Marker, sym) & strcmp(num2str(v.Color), col_rgb), H, 'uni', 0);
            G = cell2mat(G');
            for i = 1:length(list)
                t = find(G(:,i)==1);
                if ~isempty(t)
                    H{t}.Visible = 'on';
                else
                    warning([list{i} ' not found']);
                end
            end
        end

        hLegend = findall(f,'tag','legend');
        uic = get(hLegend,'UIContextMenu');
        uimenu_refresh = findall(uic,'Label','Refresh');
        callback = get(uimenu_refresh,'Callback');
        hgfeval(callback,[],[]);
end

