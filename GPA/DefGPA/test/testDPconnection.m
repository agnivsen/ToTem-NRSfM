clear all
close all
clc

numData = 5;
D = 10 * rand(3, 5);


% exact rotation
for cnt = 1 : numData
    axis = randn(3,1);
    axis = axis / norm(axis);
    Rot = SO3.Exp(0.5 * cnt * axis);
    Data{cnt} = Rot * D;
end


% whitening
for cnt = 1 : numData
   Data{cnt} = Data{cnt} + 2 * randn (size(D,1), size(D,2));
end


DT = [];
for cnt = 1 : numData
   DT = [DT, Data{cnt}'];
end


MV = [];
for cnt = 1 : numData
[~, ~, V] = svd(Data{cnt}, 'econ');
MV = [MV, V];
end


PS = zeros(size(D,2));
for cnt = 1 : numData
   PS = PS + Data{cnt}' * inv(Data{cnt} * Data{cnt}') * Data{cnt};
end
norm(PS - MV * MV')


[VDT, eigD] = byEigen (DT * DT', 3, 'largest');

[VMV, eigV] = byEigen (MV * MV', 3, 'largest');

er1 = cmpSignDiff (VDT, VMV)

eigD
eigV

[VDT, VMV]







function err = cmpSignDiff (A, B)
    
    dd = min(A-B, A+B);
    
    err = norm (dd);

end





% K, the matrix to decompose,
% d, the number of eigen vectors required.
function [V, eigValEigen, eigenValueAll] = byEigen (K, d, str)

offset = 0;

[Q, D ] = eig( (K+K')/2 );

eigVal = diag(D);

[eigVal, Index] = sort(eigVal, 'descend');  Q = Q (:, Index);

if  strcmp (str, 'largest')
    
    V = Q (:, 1+offset:d+offset);
    
    eigValEigen =  eigVal (1+offset:d+offset);
    
end

if strcmp (str, 'smallest')
    
    V = Q (:, end-d+1-offset:end-offset);
    
    eigValEigen =  eigVal (end-d+1-offset:end-offset);
    
end

eigenValueAll = eigVal;

if  norm( V' * V - eye(d) ) > 1e-12
    
    fprintf(2, 'Eigen vectors are not orthnormal!.');
    
end

end