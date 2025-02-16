function [X, Y, Z] = vector2surf(vector, n)

X = vector(1, :);
Y = vector(2, :);
Z = vector(3, :);

X = reshape(X, [], n)';
Y = reshape(Y, [], n)';
Z = reshape(Z, [], n)';