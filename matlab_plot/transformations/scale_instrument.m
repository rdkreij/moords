function [Xn, Yn, Zn] = scale_instrument(Xn, Yn, Zn, scale)

if scale == 1
    return
end
xmin = [];
xmax = [];
ymin = [];
ymax = [];
zmin = [];
zmax = [];
for ii = 1:numel(Xn)
    Xn{ii} = Xn{ii}*scale;
    Yn{ii} = Yn{ii}*scale;
    Zn{ii} = Zn{ii}*scale;
end

end