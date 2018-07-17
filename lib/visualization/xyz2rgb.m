function rgb = xyz2rgb(xyz)
% rgb = xyz2rgb(xyz)
% Xiaoyan, 2018

rgb = xyz - min(xyz, [], 1);
rgb = bsxfun(@rdivide, rgb, max(rgb, [], 1));
end