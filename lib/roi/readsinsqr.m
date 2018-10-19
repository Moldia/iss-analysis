function insqr = readsinsqr(pos, sqrlim)
% insqr = readsinsqr(pos, sqrlim)
% find reads in defined square region
% sqrlim: [xmin ymin xmax ymax]
% take reads whose position >=min and <max
% Xiaoyan, 2017


insqr = pos(:,1)>=sqrlim(1) & pos(:,2)>=sqrlim(2) &...
    pos(:,1)<sqrlim(3) & pos(:,2)<sqrlim(4);

end
