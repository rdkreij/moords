function vector = surf2vector(X, Y, Z)

X = reshape(X', 1, []);
Y = reshape(Y', 1, []);
Z = reshape(Z', 1, []);

vector = [X; Y; Z];

