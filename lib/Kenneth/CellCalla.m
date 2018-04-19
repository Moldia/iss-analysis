function h = CellCalla(MyMeans, Counts, NickNames, Colors, Zoom, Pos)

OldFig = gcf;
r = 2;
nK = size(MyMeans, 2);

sr = 0.1:0.01:.3;
ns = length(sr);
L = zeros(ns, nK);
for i=1:ns;
    s = sr(i);
    p = MyMeans*s ./ (MyMeans*s + r);
    L(i,:) = Counts*log(p) + r*sum(log(1-p));
end

[~, BestPoint] = min(abs(sr-0.2));
Probs = LogLtoP(L(BestPoint,:)');
WorthShowing = find(Probs>.05);

if nargin<5 % interactive mode if no zoom/pos argument
    %figure(4490786);
    h = pie(Probs(WorthShowing), NickNames(WorthShowing));

    for i=1:length(h)/2
        set(h(i*2-1), 'FaceColor', Colors{WorthShowing(i)});
    end
    
    figure(3490786);
    imagesc(sr, 1:nK, LogLtoP(L'));
    set(gca, 'ytick', 1:nK);
    set(gca, 'yticklabel', NickNames);
    colorbar;

    figure(OldFig)
end

if nargin>=5
    hold on
    h = pie(Probs(WorthShowing), repmat({''}, sum(WorthShowing>0)));

    for i=1:length(h)/2
        hno = (i*2-1);
%         Index = find(strcmp(h(i*2).String, NickNames), 1);
        set(h(hno), 'FaceColor', Colors{WorthShowing(i)});
        set(h(hno), 'Xdata', get(h(hno), 'Xdata')*Zoom + Pos(1));
        set(h(hno), 'Ydata', get(h(hno), 'Ydata')*Zoom + Pos(2));
        set(h(hno), 'EdgeAlpha', 0);
    end
end
    