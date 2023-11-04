clear all
close all
clc


D = 10 * rand(3, 6);


[~,~, V1] = svd(D + 0.1 * randn(3, 6), 'econ');
[~,~, V2] = svd(D + 0.1 * randn(3, 6), 'econ');
[~,~, V3] = svd(D + 0.1 * randn(3, 6), 'econ');
[~,~, V4] = svd(D + 0.1 * randn(3, 6), 'econ');


MV = [V1, V2, V3, V4]

r = rank(MV)


[U, S, ~] = svd (MV);

U = U(:, 1:3)
S = S(1:3, 1:3)

