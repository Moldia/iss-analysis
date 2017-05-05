T = readtable('SpotCoordinates_ROI.csv');

TopDist = 60;
BinSize = 1;

name = T.Name;
[uName] = unique(name);
n = zeros(length(uName),1);
n10 = n;
n20 = n;

for i=1:length(uName)
    Gene = uName{i};
    MyDots = find(strcmp(Gene,name));
    x = T.posX(MyDots);
    y = T.posY(MyDots);
    
    n(i) = length(MyDots);
    
    figure(1);
    plot(x,y, '.');
    title(Gene);

    [xSort, xOrd] = sort(x);
    [ySort, yOrd] = sort(y);
    [~,~,xPairs0] = CCG(xSort,1,TopDist,0);
    xPairs = xOrd(xPairs0);
    
    [~,~,yPairs0] = CCG(ySort,1,TopDist,0);
    yPairs = yOrd(yPairs0);
    Pairs = intersect(xPairs, yPairs, 'rows');
    
    figure(2); clf; hold on
    dx = x(Pairs(:,2))-x(Pairs(:,1));
    dy = y(Pairs(:,2))-y(Pairs(:,1));
    plot(dx, dy, '.');
    rectangle('Position', [-1 -1 2 2]*20, 'curvature', [1 1], 'edgecolor', 'r');
    title(Gene);
    
    A = accumarray(TopDist/2/BinSize+1+[ceil(dx/BinSize), ceil(dy/BinSize)], 1);
    figure(3);
    imagesc(A');
    set(gca, 'ydir', 'normal');
    title(Gene);
        rectangle('Position', [TopDist/2-20 TopDist/2-20 40 40], 'curvature', [1 1], 'edgecolor', 'r');


%     Pairs = zeros(PairSize0,2);
%     nPairs = 0;
%     for bx = 1:xBlox
%         for by = 1:yBlox
%             MyPoints = find(x >= (bx-1)*Block & x < bx*Block & y >=(by-1)*Block & y<by*Block);
%             dx = bsxfun(@minus,x(MyPoints), x(MyPoints)');
%             dy = bsxfun(@minus,y(MyPoints), y(MyPoints)');
%             if length(MyPoints)>2
%                 keyboard
%             end
%         end
%     end

    n10(i) = sum(dx.^2 + dy.^2 < 10.^2)/2/n(i);
    n20(i) = sum(dx.^2 + dy.^2 < 20.^2)/2/n(i);
   
    fprintf('%s: n %d, n10 %f, n20 %f\n', Gene, n(i), n10(i), n20(i));
    pause

end

figure(4); clf; 
loglog(n, n20, 'k.');
text(n, n20, uName);
xlabel('number of detections');
ylabel('clustering coefficient');
%axis([0 max(n) 0 max(n20)]);
