function [Xn, Yn, Zn] = rotate_instrument(Xn, Yn, Zn, ho, po, ro)

rotmat = rotation_matrix(ho, po, ro);
for ii = 1:numel(Xn)
    n = size(Xn{ii}, 1);
    vector = surf2vector(Xn{ii}, Yn{ii}, Zn{ii});
    vector = rotmat*vector;
    [Xn{ii}, Yn{ii}, Zn{ii}] = vector2surf(vector, n);
end

end