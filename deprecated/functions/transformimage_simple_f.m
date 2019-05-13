function transformimage_simple_f(original,ang,yup,xleft,scale,name,varargin)

original = imrotate(original,ang);
if yup>=0
    original = original(yup*scale+1:end,:);
else
    original = [zeros(-yup*scale,size(original,2));original];
end
if xleft>=0
    original = original(:,xleft*scale+1:end);
else
    original = [zeros(size(original,1),-xleft*scale),original];
end

if ~isempty(varargin)
    size_delta = varargin{1}-size(original);
    if size_delta(2)>=0
        original = [original,zeros(size(original,1),size_delta(2))];
    else
        original = original(:,1:varargin{1}(2));
    end
    if size_delta(1)>=0
        original = [original;zeros(size_delta(1),size(original,2))];
    else
        original = original(1:varargin{1}(1),:);
    end
end
    

imwrite(original,name,'tiff','compression','none');

end