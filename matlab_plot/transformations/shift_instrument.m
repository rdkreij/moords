function [Xn, Yn, Zn] = shift_instrument(Xn, Yn, Zn, xo, yo, zo)

for ii = 1:numel(Xn)
    Xn{ii} = Xn{ii} + xo;
    Yn{ii} = Yn{ii} + yo;
    Zn{ii} = Zn{ii} + zo;
end

end